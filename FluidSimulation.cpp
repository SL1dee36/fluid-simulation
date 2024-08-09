#include <SDL.h>
#include <SDL_opengl.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <algorithm>
#include <thread>
#include <mutex>
#include <immintrin.h> // Для SIMD инструкций

// Настройки 
const int WIDTH = 770;
const int HEIGHT = 500;
const int FPS = 5000;
const int PARTICLE_RADIUS = 7;
float GRAVITY = 0.0981f;
float REST_DENSITY = 10.0f;
float GAS_CONSTANT = 2.0f;
float VISCOSITY = 0.9f;
int PARTICLE_CREATION_RATE = 10;
float COHESION_STRENGTH = 0.2f;
float DAMPING = 0.99f;
float DRAG_COEFFICIENT = 0.1f;
const int GRAB_RADIUS = 50;
float SPRING_CONSTANT = 0.2f;
const int MAX_PARTICLES = 9000; // Ограничение на количество частиц
float ROTATION_SPEED = 20.0f; // Скорость вращения частиц
int IMPULSE_STRENGTH = 50; // Сила импульса (теперь переменная)

#define DEG2RAD 0.017453292519943295f // PI / 180


// Препятствие
const int OBSTACLE_X = 200;
const int OBSTACLE_Y = 200;
const int OBSTACLE_RADIUS = 0;

// Стены
const int WALL_THICKNESS = 0;
const SDL_Color WALL_COLOR = { 255, 255, 255, 255 };
const int PAD_X = 240;
const int PAD_Y = 20;

// Градиент цветов
const int GRADIENT_STEPS = 4096;
SDL_Color color_gradient[GRADIENT_STEPS];
const float max_speed = 16.0f;

// Количество потоков
const int NUM_THREADS = std::thread::hardware_concurrency();

// Мьютекс для синхронизации доступа к particles
std::mutex particles_mutex;

// Флаг для отображения стрелок
bool show_arrows = false;

// Переменная для хранения состояния дебаг-меню
bool debug_mode = false;

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

        // Добавлено условие для верхней границы
        if (y < PARTICLE_RADIUS + PAD_Y) {
            y = PARTICLE_RADIUS + PAD_Y;
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

    void draw() {
        // Использование градиента для визуализации скорости
        float speed = std::sqrt(vx * vx + vy * vy);
        int colorIndex = static_cast<int>(speed / max_speed * (GRADIENT_STEPS - 1));
        colorIndex = std::min(colorIndex, GRADIENT_STEPS - 1);

        // Рисование круга с помощью OpenGL
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glColor4ub(color_gradient[colorIndex].r, color_gradient[colorIndex].g, color_gradient[colorIndex].b, 128);
        glBegin(GL_POLYGON);

        // Выбор количества сегментов в зависимости от debug_mode
        int num_segments = debug_mode ? 0 : 32;

        for (int i = 0; i < num_segments; i++) {
            float degInRad = i * 2 * M_PI / num_segments;
            glVertex2f(x + cos(degInRad) * PARTICLE_RADIUS, y + sin(degInRad) * PARTICLE_RADIUS);
        }
        glEnd();

        glDisable(GL_BLEND);
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

    // Применение импульса к частице
    void apply_impulse(float impulse_x, float impulse_y) {
        vx += impulse_x / density;
        vy += impulse_y / density;
    }
};

// Оптимизированная функция для нахождения соседей
// Используем SIMD для ускорения вычислений расстояний
void find_neighbors(std::vector<Particle>& particles, int start, int end) {
    for (size_t i = start; i < end; ++i) {
        particles[i].near_particles.clear();
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i == j) continue;

            // SIMD оптимизация: вычисляем 4 расстояния одновременно
            __m128 px = _mm_set1_ps(particles[i].x);
            __m128 py = _mm_set1_ps(particles[i].y);
            __m128 ox = _mm_set1_ps(particles[j].x);
            __m128 oy = _mm_set1_ps(particles[j].y);

            __m128 dx = _mm_sub_ps(px, ox);
            __m128 dy = _mm_sub_ps(py, oy);

            __m128 distance_sq = _mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy));

            // Сравниваем с квадратом суммы радиусов
            __m128 radius_sum_sq = _mm_set1_ps((2 * PARTICLE_RADIUS) * (2 * PARTICLE_RADIUS));
            __m128 cmp = _mm_cmplt_ps(distance_sq, radius_sum_sq);

            // Если хоть один элемент сравнения истинный, добавляем в соседей
            int mask = _mm_movemask_ps(cmp);
            if (mask) {
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

// Структура для хранения данных ползунка
struct Slider {
    std::string label;
    float* value;
    float min_value;
    float max_value;
    int x, y, width, height;
};

// Функция для отрисовки ползунка
void draw_slider(Slider slider) {
    // Отрисовка фона ползунка
    glColor4ub(200, 200, 200, 255);
    glBegin(GL_QUADS);
    glVertex2f(slider.x, slider.y);
    glVertex2f(slider.x + slider.width, slider.y);
    glVertex2f(slider.x + slider.width, slider.y + slider.height);
    glVertex2f(slider.x, slider.y + slider.height);
    glEnd();

    // Вычисление позиции ползунка
    float slider_pos = slider.x + slider.width * (*slider.value - slider.min_value) / (slider.max_value - slider.min_value);

    // Отрисовка ползунка
    glColor4ub(100, 100, 100, 255);
    glBegin(GL_QUADS);
    glVertex2f(slider_pos - 5, slider.y);
    glVertex2f(slider_pos + 5, slider.y);
    glVertex2f(slider_pos + 5, slider.y + slider.height);
    glVertex2f(slider_pos - 5, slider.y + slider.height);
    glEnd();

    // Отрисовка текста
    // (Здесь вам нужно использовать библиотеку для отрисовки текста, 
    //  например, SDL_ttf. Этот код показывает только позицию текста)
    // int text_x = slider.x + slider.width / 2;
    // int text_y = slider.y + slider.height + 10;
    // render_text(slider.label, text_x, text_y);
}

// Функция для обработки событий ползунка
bool handle_slider_event(Slider slider, int mouseX, int mouseY) {
    if (mouseX >= slider.x && mouseX <= slider.x + slider.width &&
        mouseY >= slider.y && mouseY <= slider.y + slider.height) {
        *slider.value = slider.min_value + (mouseX - slider.x) * (slider.max_value - slider.min_value) / slider.width;
        return true;
    }
    return false;
}

int main(int argc, char* argv[]) {
    // Инициализация SDL
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cerr << "SDL_Init Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    // Настройка OpenGL
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 32);

    // Создание окна
    SDL_Window* window = SDL_CreateWindow("2D Fluid Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        WIDTH, HEIGHT, SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL);
    if (window == nullptr) {
        std::cerr << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return 1;
    }

    // Создание контекста OpenGL
    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    if (gl_context == nullptr) {
        std::cerr << "SDL_GL_CreateContext Error: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    // Инициализация OpenGL
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, WIDTH, HEIGHT, 0, 1, -1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Создание градиента
    SDL_Color color1 = { 0, 0, 255, 255 }; // Синий
    SDL_Color color2 = { 255, 255, 255, 255 }; // Желтый
    create_gradient(color1, color2, GRADIENT_STEPS, color_gradient);

    // Создание списка частиц
    std::vector<Particle> particles;

    // Флаги для кнопок мыши
    bool mouse_pressed_left = false;
    bool mouse_pressed_right = false;
    bool mouse_pressed_central = false;

    // Список захваченных частиц
    std::vector<Particle*> grabbed_particles;

    // Генератор случайных чисел
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(-PARTICLE_RADIUS * 2, PARTICLE_RADIUS * 2);

    // Создание ползунков
    std::vector<Slider> sliders = {
        {"Gravity", &GRAVITY, 0.0f, 1.0f, 10, 10, 100, 20},
        {"Rest Density", &REST_DENSITY, 1.0f, 20.0f, 10, 40, 100, 20},
        {"Gas Constant", &GAS_CONSTANT, 0.1f, 5.0f, 10, 70, 100, 20},
        {"Viscosity", &VISCOSITY, 0.1f, 2.0f, 10, 100, 100, 20},
        {"Cohesion", &COHESION_STRENGTH, 0.0f, 1.0f, 10, 130, 100, 20},
        {"Spring Constant", &SPRING_CONSTANT, 0.0f, 1.0f, 10, 160, 100, 20},
        // Добавьте остальные ползунки здесь...
    };

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

                    // Обработка событий ползунков
                    int mouseX, mouseY;
                    SDL_GetMouseState(&mouseX, &mouseY);
                    for (Slider& slider : sliders) {
                        handle_slider_event(slider, mouseX, mouseY);
                    }
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
                // Обработка нажатия средней кнопки мыши
                else if (event.button.button == SDL_BUTTON_MIDDLE) {
                    mouse_pressed_central = true;
                    int mouseX, mouseY;
                    SDL_GetMouseState(&mouseX, &mouseY);
                    // Применение импульса к частицам в радиусе действия
                    for (Particle& particle : particles) {
                        float dx = particle.x - mouseX;
                        float dy = particle.y - mouseY;
                        float distance_sq = dx * dx + dy * dy;
                        float impulse_radius_sq = GRAB_RADIUS * GRAB_RADIUS;

                        if (distance_sq <= impulse_radius_sq) {
                            float distance = std::sqrt(distance_sq);
                            // Расчет импульса в зависимости от расстояния
                            float impulse_strength = IMPULSE_STRENGTH * (1 - distance / GRAB_RADIUS);
                            float impulse_x = impulse_strength * dx / distance;
                            float impulse_y = impulse_strength * dy / distance;
                            particle.apply_impulse(impulse_x, impulse_y);
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
                else if (event.button.button == SDL_BUTTON_MIDDLE) {
                    mouse_pressed_central = false;
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
                // Управление силой импульса с помощью клавиш
                else if (event.key.keysym.sym == SDLK_UP) {
                    IMPULSE_STRENGTH += 5.0f;
                }
                else if (event.key.keysym.sym == SDLK_DOWN) {
                    IMPULSE_STRENGTH = std::max(0.0f, IMPULSE_STRENGTH - 5.0f);
                }
                // Включение/выключение стрелок при нажатии на 1 - Дебаг меню
                else if (event.key.keysym.sym == SDLK_1) {
                    debug_mode = !debug_mode;
                    show_arrows = !show_arrows;
                }
                // При нажатии 0 - const int NUM_SEGMENTS = 32
                else if (event.key.keysym.sym == SDLK_0) {
                    debug_mode = false; // Выключаем дебаг-режим
                    show_arrows = false;
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

        if (mouse_pressed_central) {
            int mouseX, mouseY;
            SDL_GetMouseState(&mouseX, &mouseY);
            // Применение импульса к частицам в радиусе действия
            for (Particle& particle : particles) {
                float dx = particle.x - mouseX;
                float dy = particle.y - mouseY;
                float distance_sq = dx * dx + dy * dy;
                float impulse_radius_sq = GRAB_RADIUS * GRAB_RADIUS;

                if (distance_sq <= impulse_radius_sq) {
                    float distance = std::sqrt(distance_sq);
                    // Расчет импульса в зависимости от расстояния
                    float impulse_strength = IMPULSE_STRENGTH * (1 - distance / GRAB_RADIUS);
                    float impulse_x = impulse_strength * dx / distance;
                    float impulse_y = impulse_strength * dy / distance;
                    particle.apply_impulse(impulse_x, impulse_y);
                }
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
        glClear(GL_COLOR_BUFFER_BIT);

        // Рисуем стены
        glColor4ub(WALL_COLOR.r, WALL_COLOR.g, WALL_COLOR.b, WALL_COLOR.a);
        glBegin(GL_LINES);
        glVertex2f(PAD_X, PAD_Y);
        glVertex2f(WIDTH - PAD_X, PAD_Y);
        glVertex2f(PAD_X, HEIGHT - PAD_Y);
        glVertex2f(WIDTH - PAD_X, HEIGHT - PAD_Y);
        glVertex2f(PAD_X, PAD_Y);
        glVertex2f(PAD_X, HEIGHT - PAD_Y);
        glVertex2f(WIDTH - PAD_X, PAD_Y);
        glVertex2f(WIDTH - PAD_X, HEIGHT - PAD_Y);
        glEnd();

        // Рисуем препятствие
        glColor4ub(255, 255, 255, 255);
        glBegin(GL_POLYGON);
        const int NUM_SEGMENTS_OBSTACLE = 365; // Уменьшаем количество сегментов для препятствия
        for (int i = 0; i < NUM_SEGMENTS_OBSTACLE; i++) {
            float degInRad = i * 2 * M_PI / NUM_SEGMENTS_OBSTACLE;
            glVertex2f(OBSTACLE_X + cos(degInRad) * OBSTACLE_RADIUS, OBSTACLE_Y + sin(degInRad) * OBSTACLE_RADIUS);
        }
        glEnd();

        particles_mutex.lock();
        for (Particle& particle : particles) {
            particle.draw();
        }
        particles_mutex.unlock();

        // Рисуем стрелки поверх частиц
        particles_mutex.lock();
        for (Particle& particle : particles) {
            if (show_arrows) {
                float speed = std::sqrt(particle.vx * particle.vx + particle.vy * particle.vy);
                if (speed > 0.1f) { // Рисуем стрелку только если скорость значительна
                    float arrow_length = 10.0f;
                    float arrow_angle = atan2(particle.vy, particle.vx);
                    float arrow_tip_x = particle.x + arrow_length * cos(arrow_angle);
                    float arrow_tip_y = particle.y + arrow_length * sin(arrow_angle);
                    glColor4ub(255, 0, 0, 255); // Красный цвет для стрелки
                    glBegin(GL_LINES);
                    glVertex2f(particle.x, particle.y);
                    glVertex2f(arrow_tip_x, arrow_tip_y);
                    glEnd();
                }
            }
        }
        particles_mutex.unlock();

        // Отрисовка ползунков
        for (Slider& slider : sliders) {
            draw_slider(slider);
        }

        // Вывод FPS и количества частиц
        Uint32 current_time = SDL_GetTicks();
        float fps = 1000.0f / (current_time - last_time);
        last_time = current_time;

        // Преобразование FPS в целое число
        int fps_int = static_cast<int>(fps);

        std::string title = "2D Fluid Simulation v0.975R | Particles: " + std::to_string(particles.size()) + "/" + std::to_string(MAX_PARTICLES) + " | R - clean up | FPS: " + std::to_string(fps_int) + " | Impulse Strength: " + std::to_string(IMPULSE_STRENGTH);
        SDL_SetWindowTitle(window, title.c_str());

        SDL_GL_SwapWindow(window);
        SDL_Delay(1000 / FPS);
    }

    // Освобождение ресурсов
    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}