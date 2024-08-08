#include <SDL.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <algorithm>
#include <thread>
#include <mutex>

// Настройки 
const int WIDTH = 600;
const int HEIGHT = 720;
const int FPS = 2400;
const int PARTICLE_RADIUS = 3;
const float GRAVITY = 0.0981f;
const float REST_DENSITY = 0.00001f;
const float GAS_CONSTANT = 2.0;
const float VISCOSITY = 0.9f;
const int PARTICLE_CREATION_RATE = 10;
const float COHESION_STRENGTH = 0.1f;
const float DAMPING = 0.99f;
const float DRAG_COEFFICIENT = 0.1f;
const int GRAB_RADIUS = 100;
const float SPRING_CONSTANT = 0.0f;
const int MAX_PARTICLES = 5000; // Ограничение на количество частиц
const float ROTATION_SPEED = 0.05f; // Скорость вращения частиц

// Препятствие
const int OBSTACLE_X = 400;
const int OBSTACLE_Y = 300;
const int OBSTACLE_RADIUS = 100;

// Стены
const int WALL_THICKNESS = 0;
const SDL_Color WALL_COLOR = { 255, 255, 255, 255 };
const int PAD_X = 20;
const int PAD_Y = 20;

// Градиент цветов
const int GRADIENT_STEPS = 256;
SDL_Color color_gradient[GRADIENT_STEPS];
const float max_speed = 16.0f;

// Количество потоков
const int NUM_THREADS = std::thread::hardware_concurrency();

// Мьютекс для синхронизации доступа к particles
std::mutex particles_mutex;

// Функция для создания градиента цветов
void create_gradient(SDL_Color color1, SDL_Color color2, int steps, SDL_Color* gradient) {
    for (int i = 0; i < steps; ++i) {
        float t = i / (steps - 1.0f);
        gradient[i].r = static_cast<Uint8>(color1.r * (1 - t) + color2.r * t);
        gradient[i].g = static_cast<Uint8>(color1.g * (1 - t) + color2.g * t);
        gradient[i].b = static_cast<Uint8>(color1.b * (1 - t) + color2.b * t);
        gradient[i].a = 255;
    }
}

// Функция для рисования заполненного круга
void filledCircleRGBA(SDL_Renderer* renderer, int x, int y, int radius, Uint8 r, Uint8 g, Uint8 b, Uint8 a) {
    for (int dy = -radius; dy <= radius; ++dy) {
        for (int dx = -radius; dx <= radius; ++dx) {
            if (dx * dx + dy * dy <= radius * radius) {
                SDL_SetRenderDrawColor(renderer, r, g, b, a);
                SDL_RenderDrawPoint(renderer, x + dx, y + dy);
            }
        }
    }
}

// Класс частицы
class Particle {
public:
    float x;
    float y;
    float vx;
    float vy;
    float density;
    float pressure;
    float angle; // Угол поворота
    float angular_velocity; // Угловая скорость
    std::vector<Particle*> near_particles;
    bool grabbed;

    Particle(float x, float y) : x(x), y(y), vx(0), vy(0), density(0), pressure(0), angle(0), angular_velocity(0), grabbed(false) {}

    void update() {
        vy += GRAVITY;
        x += vx;
        y += vy;

        // Обновление вращения
        angle += angular_velocity;

        // Столкновения с границами (оптимизировано)
        if (x < PARTICLE_RADIUS + PAD_X) {
            x = PARTICLE_RADIUS + PAD_X;
            vx *= -0.7f;
        }
        else if (x > WIDTH - PARTICLE_RADIUS - PAD_X) {
            x = WIDTH - PARTICLE_RADIUS - PAD_X;
            vx *= -0.7f;
        }

        if (y > HEIGHT - PARTICLE_RADIUS - PAD_Y) {
            y = HEIGHT - PARTICLE_RADIUS - PAD_Y;
            vy *= -0.7f;
        }

        // Столкновение с препятствием (оптимизировано)
        float dx = x - OBSTACLE_X;
        float dy = y - OBSTACLE_Y;
        float distance_to_obstacle_sq = dx * dx + dy * dy;
        float radius_sum_sq = (PARTICLE_RADIUS + OBSTACLE_RADIUS) * (PARTICLE_RADIUS + OBSTACLE_RADIUS);

        if (distance_to_obstacle_sq < radius_sum_sq) {
            float distance_to_obstacle = std::sqrt(distance_to_obstacle_sq);
            float overlap = (PARTICLE_RADIUS + OBSTACLE_RADIUS) - distance_to_obstacle;
            x += overlap * dx / distance_to_obstacle;
            y += overlap * dy / distance_to_obstacle;
            vx *= -0.7f;
            vy *= -0.7f;
        }

        // Демпфирование скорости
        vx *= DAMPING;
        vy *= DAMPING;
    }

    void draw(SDL_Renderer* renderer) {
        // Использование градиента
        float speed = std::sqrt(vx * vx + vy * vy);
        int colorIndex = static_cast<int>(speed / max_speed * (GRADIENT_STEPS - 1));
        colorIndex = std::min(colorIndex, GRADIENT_STEPS - 1);

        // Рисование круга с учетом вращения
        // (Здесь можно добавить более сложную логику для визуализации вращения)
        filledCircleRGBA(renderer, static_cast<int>(x), static_cast<int>(y), PARTICLE_RADIUS,
            color_gradient[colorIndex].r, color_gradient[colorIndex].g,
            color_gradient[colorIndex].b, color_gradient[colorIndex].a);
    }

    void reset_velocity() {
        vx = 0;
        vy = 0;
        angular_velocity = 0; // Сброс угловой скорости
    }

    void apply_spring_force(Particle* other) {
        float dx = x - other->x;
        float dy = y - other->y;
        float distance_sq = dx * dx + dy * dy;
        float radius_sum_sq = (2 * PARTICLE_RADIUS) * (2 * PARTICLE_RADIUS);

        if (distance_sq > radius_sum_sq) {
            float distance = std::sqrt(distance_sq);
            float force = SPRING_CONSTANT * (distance - 2 * PARTICLE_RADIUS);
            float fx = force * dx / distance;
            float fy = force * dy / distance;
            vx -= fx / density;
            vy -= fy / density;
            other->vx += fx / other->density;
            other->vy += fy / other->density;

            // Добавление вращения
            float torque = (fx * dy - fy * dx) * ROTATION_SPEED;
            angular_velocity += torque / density;
            other->angular_velocity -= torque / other->density;
        }
    }
};

// Оптимизированная функция для нахождения соседей
void find_neighbors(std::vector<Particle>& particles, int start, int end) {
    for (size_t i = start; i < end; ++i) {
        particles[i].near_particles.clear();
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i == j) continue;
            float dx = particles[i].x - particles[j].x;
            float dy = particles[i].y - particles[j].y;
            float distance_sq = dx * dx + dy * dy;
            float radius_sum_sq = (2 * PARTICLE_RADIUS) * (2 * PARTICLE_RADIUS);

            if (distance_sq < radius_sum_sq) {
                particles[i].near_particles.push_back(&particles[j]);
            }
        }
    }
}

// Оптимизированная функция для вычисления плотности и давления
void calculate_density_pressure(std::vector<Particle>& particles, int start, int end) {
    for (size_t i = start; i < end; ++i) {
        particles[i].density = 0;
        for (Particle* other : particles[i].near_particles) {
            float dx = particles[i].x - other->x;
            float dy = particles[i].y - other->y;
            float distance_sq = dx * dx + dy * dy;
            float radius_sum_sq = (2 * PARTICLE_RADIUS) * (2 * PARTICLE_RADIUS);

            if (distance_sq < radius_sum_sq) {
                particles[i].density += 1;
            }
        }
        particles[i].density = std::max(particles[i].density, 0.1f);
        particles[i].pressure = GAS_CONSTANT * (particles[i].density - REST_DENSITY);
    }
}

// Оптимизированная функция для вычисления сил
void calculate_forces(std::vector<Particle>& particles, int start, int end) {
    for (size_t i = start; i < end; ++i) {
        float dx = 0, dy = 0;
        for (Particle* other : particles[i].near_particles) {
            float dx_diff = particles[i].x - other->x;
            float dy_diff = particles[i].y - other->y;
            float distance_sq = dx_diff * dx_diff + dy_diff * dy_diff;
            float radius_sum_sq = (2 * PARTICLE_RADIUS) * (2 * PARTICLE_RADIUS);

            if (distance_sq < radius_sum_sq && distance_sq > 0) {
                float distance = std::sqrt(distance_sq);

                // Сила давления
                float pressure_force = (particles[i].pressure + other->pressure) / 2 * (1 - distance / (2 * PARTICLE_RADIUS));
                dx += pressure_force * dx_diff / distance;
                dy += pressure_force * dy_diff / distance;

                // Сила вязкости
                float vx_diff = other->vx - particles[i].vx;
                float vy_diff = other->vy - particles[i].vy;
                float viscosity_force = VISCOSITY * (vx_diff * dx_diff + vy_diff * dy_diff) / distance_sq;
                dx += viscosity_force * dx_diff;
                dy += viscosity_force * dy_diff;

                // Сила сцепления
                float cohesion_force = COHESION_STRENGTH / distance;
                dx -= cohesion_force * dx_diff;
                dy -= cohesion_force * dy_diff;

                // Сила пружины
                particles[i].apply_spring_force(other);
            }
        }

        // Применение сил
        particles[i].vx += dx / particles[i].density;
        particles[i].vy += dy / particles[i].density;
    }
}

// Функция для обновления физики в отдельном потоке
void update_physics_thread(std::vector<Particle>& particles, int start, int end) {
    find_neighbors(particles, start, end);
    calculate_density_pressure(particles, start, end);
    calculate_forces(particles, start, end);
}

int main(int argc, char* argv[]) {
    // Инициализация SDL
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cerr << "SDL_Init Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    // Создание окна
    SDL_Window* window = SDL_CreateWindow("2D Fluid Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
    if (window == nullptr) {
        std::cerr << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return 1;
    }

    // Создание рендерера
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == nullptr) {
        std::cerr << "SDL_CreateRenderer Error: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    // Создание градиента
    SDL_Color color1 = { 0, 0, 255, 255 };
    SDL_Color color2 = { 255, 0, 0, 255 };
    create_gradient(color1, color2, GRADIENT_STEPS, color_gradient);

    // Создание списка частиц
    std::vector<Particle> particles;

    // Флаги для кнопок мыши
    bool mouse_pressed_left = false;
    bool mouse_pressed_right = false;

    // Список захваченных частиц
    std::vector<Particle*> grabbed_particles;

    // Генератор случайных чисел
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(-PARTICLE_RADIUS * 2, PARTICLE_RADIUS * 2);

    // Главный цикл
    bool running = true;
    Uint32 last_time = SDL_GetTicks();
    while (running) {
        // Обработка событий
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
            else if (event.type == SDL_MOUSEBUTTONDOWN) {
                if (event.button.button == SDL_BUTTON_LEFT) {
                    mouse_pressed_left = true;
                }

                else if (event.button.button == SDL_BUTTON_RIGHT) {
                    mouse_pressed_right = true;
                    int mouseX, mouseY;
                    SDL_GetMouseState(&mouseX, &mouseY);
                    for (Particle& particle : particles) {
                        float dx = particle.x - mouseX;
                        float dy = particle.y - mouseY;
                        float distance_sq = dx * dx + dy * dy;
                        float grab_radius_sq = GRAB_RADIUS * GRAB_RADIUS;

                        if (distance_sq <= grab_radius_sq) {
                            particle.grabbed = true;
                            grabbed_particles.push_back(&particle);
                        }
                    }
                }
            }
            else if (event.type == SDL_MOUSEBUTTONUP) {
                if (event.button.button == SDL_BUTTON_LEFT) {
                    mouse_pressed_left = false;
                }
                else if (event.button.button == SDL_BUTTON_RIGHT) {
                    mouse_pressed_right = false;
                    for (Particle* particle : grabbed_particles) {
                        particle->grabbed = false;
                    }
                    grabbed_particles.clear();
                }
            }
            else if (event.type == SDL_KEYDOWN) {
                if (event.key.keysym.sym == SDLK_s) {
                    for (Particle& particle : particles) {
                        particle.reset_velocity();
                    }
                }
                else if (event.key.keysym.sym == SDLK_r) { // Удаление всех частиц
                    particles.clear();
                }
            }
        }

        // Создание новых частиц
        if (mouse_pressed_left && particles.size() < MAX_PARTICLES) {
            for (int i = 0; i < PARTICLE_CREATION_RATE; ++i) {
                int mouseX, mouseY;
                SDL_GetMouseState(&mouseX, &mouseY);
                float x = mouseX + distrib(gen);
                float y = mouseY + distrib(gen);
                particles_mutex.lock();
                particles.emplace_back(x, y);
                particles_mutex.unlock();
            }
        }

        // Перемещение захваченных частиц
        if (mouse_pressed_right) {
            int mouseX, mouseY;
            SDL_GetMouseState(&mouseX, &mouseY);
            for (Particle* particle : grabbed_particles) {
                float dx = mouseX - particle->x;
                float dy = mouseY - particle->y;
                particle->vx = dx * DRAG_COEFFICIENT;
                particle->vy = dy * DRAG_COEFFICIENT;
            }
        }

        // Обновление физики в нескольких потоках
        std::vector<std::thread> threads;
        int particles_per_thread = particles.size() / NUM_THREADS;
        for (int i = 0; i < NUM_THREADS; ++i) {
            int start = i * particles_per_thread;
            int end = (i == NUM_THREADS - 1) ? particles.size() : (i + 1) * particles_per_thread;
            threads.push_back(std::thread(update_physics_thread, std::ref(particles), start, end));
        }

        // Ожидание завершения всех потоков
        for (auto& thread : threads) {
            thread.join();
        }

        // Обновление частиц
        particles_mutex.lock();
        for (Particle& particle : particles) {
            particle.update();
        }
        particles_mutex.unlock();

        // Отрисовка
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // Рисуем стены (оптимизировано - толщина 0)
        SDL_SetRenderDrawColor(renderer, WALL_COLOR.r, WALL_COLOR.g, WALL_COLOR.b, WALL_COLOR.a);
        SDL_RenderDrawLine(renderer, PAD_X, PAD_Y, WIDTH - PAD_X, PAD_Y);
        SDL_RenderDrawLine(renderer, PAD_X, HEIGHT - PAD_Y, WIDTH - PAD_X, HEIGHT - PAD_Y);
        SDL_RenderDrawLine(renderer, PAD_X, PAD_Y, PAD_X, HEIGHT - PAD_Y);
        SDL_RenderDrawLine(renderer, WIDTH - PAD_X, PAD_Y, WIDTH - PAD_X, HEIGHT - PAD_Y);

        // Рисуем препятствие (оптимизировано - рисование круга)
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        filledCircleRGBA(renderer, OBSTACLE_X, OBSTACLE_Y, OBSTACLE_RADIUS, 255, 255, 255, 255);

        // Рисуем частицы
        particles_mutex.lock();
        for (Particle& particle : particles) {
            particle.draw(renderer);
        }
        particles_mutex.unlock();

        // Вывод FPS и количества частиц
        Uint32 current_time = SDL_GetTicks();
        float fps = 1000.0f / (current_time - last_time);
        last_time = current_time;

        // Преобразование FPS в целое число
        int fps_int = static_cast<int>(fps);

        std::string title = "2D Fluid Simulation v0.910 | Particles: " + std::to_string(particles.size()) + "/" + std::to_string(MAX_PARTICLES) + " | R - clean up | FPS: " + std::to_string(fps_int);
        SDL_SetWindowTitle(window, title.c_str());

        SDL_RenderPresent(renderer);
        SDL_Delay(1000 / FPS);
    }

    // Освобождение ресурсов
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}