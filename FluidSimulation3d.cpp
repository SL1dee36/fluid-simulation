
#include <SDL.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <algorithm>

// ���������
const int WIDTH = 600;
const int HEIGHT = 720;
const int FPS = 6000;
const float PARTICLE_RADIUS = 10.0f;
const float GRAVITY = 0.0981f;
const float REST_DENSITY = 0.01f;
const float GAS_CONSTANT = 0.95f;
const float VISCOSITY = 0.9f;
const int PARTICLE_CREATION_RATE = 10;
const float COHESION_STRENGTH = 0.0f;
const float DAMPING = 0.99f;
const float DRAG_COEFFICIENT = 0.1f;
const int GRAB_RADIUS = 100;
const float SPRING_CONSTANT = 0.25f;
const int MAX_PARTICLES = 150000; // ����������� �� ���������� ������

// ����������� (������ �����)
const float OBSTACLE_X = 400.0f;
const float OBSTACLE_Y = 300.0f;
const float OBSTACLE_Z = 200.0f;
const float OBSTACLE_RADIUS = 100.0f;

// ����� (������ �������������� ���)
const float WALL_THICKNESS = 0.0f;
const SDL_Color WALL_COLOR = { 255, 255, 255, 255 };
const float PAD_X = 20.0f;
const float PAD_Y = 20.0f;
const float PAD_Z = 20.0f;

// �������� ������
const int GRADIENT_STEPS = 256;
SDL_Color color_gradient[GRADIENT_STEPS];
const float max_speed = 25.0f;

// ���� ������
float cameraYaw = 0.0f; // ������� ������ ��� Y
float cameraPitch = 0.0f; // ������ �����/����

// ������� ������
float cameraX = 0.0f;
float cameraY = 0.0f;
float cameraZ = 0.0f;

// ������� ��� �������� ��������� ������ (�� ����������)
void create_gradient(SDL_Color color1, SDL_Color color2, int steps, SDL_Color* gradient) {
    for (int i = 0; i < steps; ++i) {
        float t = i / (steps - 1.0f);
        gradient[i].r = static_cast<Uint8>(color1.r * (1 - t) + color2.r * t);
        gradient[i].g = static_cast<Uint8>(color1.g * (1 - t) + color2.g * t);
        gradient[i].b = static_cast<Uint8>(color1.b * (1 - t) + color2.b * t);
        gradient[i].a = 255;
    }
}

// ������� ��� ��������� ������������ ����� (������ ������ � ������ �����������, ����� ������ � ������� ������)
void filledCircleRGBA(SDL_Renderer* renderer, float x, float y, float z, float radius, Uint8 r, Uint8 g, Uint8 b, Uint8 a) {
    // ���������� ������� ������
    x -= cameraX;
    y -= cameraY;
    z -= cameraZ;

    // ������� ������ ��� Y (Yaw)
    float rotatedX = x * cos(cameraYaw) - z * sin(cameraYaw);
    float rotatedZ = x * sin(cameraYaw) + z * cos(cameraYaw);

    // ������ �����/���� (Pitch) - �� ����������, ������� ����� ������� ������� ��������

    // ������� �������� ����������� (����� ��������)
    float perspectiveScale = 1.0f / (1.0f + rotatedZ / 500.0f);
    int screenX = static_cast<int>(rotatedX * perspectiveScale + WIDTH / 2);
    int screenY = static_cast<int>(y * perspectiveScale + HEIGHT / 2);
    int screenRadius = static_cast<int>(radius * perspectiveScale);

    for (int dy = -screenRadius; dy <= screenRadius; ++dy) {
        for (int dx = -screenRadius; dx <= screenRadius; ++dx) {
            if (dx * dx + dy * dy <= screenRadius * screenRadius) {
                SDL_SetRenderDrawColor(renderer, r, g, b, a);
                SDL_RenderDrawPoint(renderer, screenX + dx, screenY + dy);
            }
        }
    }
}

// ����� ������� (�������� z)
class Particle {
public:
    float x;
    float y;
    float z; // ����� ����������
    float vx;
    float vy;
    float vz; // ����� ��������
    float density;
    float pressure;
    std::vector<Particle*> near_particles;
    bool grabbed;

    Particle(float x, float y, float z) : x(x), y(y), z(z), vx(0), vy(0), vz(0), density(0), pressure(0), grabbed(false) {}

    void update() {
        vy += GRAVITY;
        x += vx;
        y += vy;
        z += vz; // ���������� z

        // ������������ � ��������� (���)
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
        // ��������� ����������� �� z
        if (z < PARTICLE_RADIUS + PAD_Z) {
            z = PARTICLE_RADIUS + PAD_Z;
            vz *= -0.7f;
        }
        else if (z > 400 - PARTICLE_RADIUS - PAD_Z) { // �������������� ������������ �������
            z = 400 - PARTICLE_RADIUS - PAD_Z;
            vz *= -0.7f;
        }


        // ������������ � ������������ (�����)
        float dx = x - OBSTACLE_X;
        float dy = y - OBSTACLE_Y;
        float dz = z - OBSTACLE_Z;
        float distance_to_obstacle_sq = dx * dx + dy * dy + dz * dz; // 3D ����������
        float radius_sum_sq = (PARTICLE_RADIUS + OBSTACLE_RADIUS) * (PARTICLE_RADIUS + OBSTACLE_RADIUS);

        if (distance_to_obstacle_sq < radius_sum_sq) {
            float distance_to_obstacle = std::sqrt(distance_to_obstacle_sq);
            float overlap = (PARTICLE_RADIUS + OBSTACLE_RADIUS) - distance_to_obstacle;
            x += overlap * dx / distance_to_obstacle;
            y += overlap * dy / distance_to_obstacle;
            z += overlap * dz / distance_to_obstacle; // ���������� z
            vx *= -0.7f;
            vy *= -0.7f;
            vz *= -0.7f; // ���������� vz
        }

        // ������������� ��������
        vx *= DAMPING;
        vy *= DAMPING;
        vz *= DAMPING; // ������������� vz
    }

    void draw(SDL_Renderer* renderer) {
        // ������������� ��������� (�� ����������)
        float speed = std::sqrt(vx * vx + vy * vy + vz * vz); // 3D speed
        int colorIndex = static_cast<int>(speed / max_speed * (GRADIENT_STEPS - 1));
        colorIndex = std::min(colorIndex, GRADIENT_STEPS - 1);

        // ��������� ����� � ������ �����������
        filledCircleRGBA(renderer, x, y, z, PARTICLE_RADIUS,
            color_gradient[colorIndex].r, color_gradient[colorIndex].g,
            color_gradient[colorIndex].b, color_gradient[colorIndex].a);
    }

    void reset_velocity() {
        vx = 0;
        vy = 0;
        vz = 0; // ����� vz
    }

    void apply_spring_force(Particle* other) {
        float dx = x - other->x;
        float dy = y - other->y;
        float dz = z - other->z; // ������� �� z
        float distance_sq = dx * dx + dy * dy + dz * dz; // 3D ����������
        float radius_sum_sq = (2 * PARTICLE_RADIUS) * (2 * PARTICLE_RADIUS);

        if (distance_sq > radius_sum_sq) {
            float distance = std::sqrt(distance_sq);
            float force = SPRING_CONSTANT * (distance - 2 * PARTICLE_RADIUS);
            float fx = force * dx / distance;
            float fy = force * dy / distance;
            float fz = force * dz / distance; // ���� �� z
            vx -= fx / density;
            vy -= fy / density;
            vz -= fz / density; // ���������� ���� �� z
            other->vx += fx / other->density;
            other->vy += fy / other->density;
            other->vz += fz / other->density; // ���������� ���� �� z
        }
    }
};

// ���������������� ������� ��� ���������� ������� (�������� ��� 3D)
void find_neighbors(std::vector<Particle>& particles) {
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].near_particles.clear();
        for (size_t j = i + 1; j < particles.size(); ++j) {
            float dx = particles[i].x - particles[j].x;
            float dy = particles[i].y - particles[j].y;
            float dz = particles[i].z - particles[j].z; // ������� �� z
            float distance_sq = dx * dx + dy * dy + dz * dz; // 3D ����������
            float radius_sum_sq = (2 * PARTICLE_RADIUS) * (2 * PARTICLE_RADIUS);

            if (distance_sq < radius_sum_sq) {
                particles[i].near_particles.push_back(&particles[j]);
                particles[j].near_particles.push_back(&particles[i]);
            }
        }
    }
}

// ���������������� ������� ��� ���������� ��������� � �������� (�� ����������)
void calculate_density_pressure(std::vector<Particle>& particles) {
    for (Particle& particle : particles) {
        particle.density = 0;
        for (Particle* other : particle.near_particles) {
            float dx = particle.x - other->x;
            float dy = particle.y - other->y;
            float dz = particle.z - other->z; // ������� �� z (�� ������������)
            float distance_sq = dx * dx + dy * dy + dz * dz; // 3D ����������
            float radius_sum_sq = (2 * PARTICLE_RADIUS) * (2 * PARTICLE_RADIUS);

            if (distance_sq < radius_sum_sq) {
                particle.density += 1;
            }
        }
        particle.density = std::max(particle.density, 0.1f);
        particle.pressure = GAS_CONSTANT * (particle.density - REST_DENSITY);
    }
}

// ���������������� ������� ��� ���������� ��� (�������� ��� 3D)
void calculate_forces(std::vector<Particle>& particles) {
    for (Particle& particle : particles) {
        float dx = 0, dy = 0, dz = 0; // �������� dz
        for (Particle* other : particle.near_particles) {
            float dx_diff = particle.x - other->x;
            float dy_diff = particle.y - other->y;
            float dz_diff = particle.z - other->z; // ������� �� z
            float distance_sq = dx_diff * dx_diff + dy_diff * dy_diff + dz_diff * dz_diff; // 3D ����������
            float radius_sum_sq = (2 * PARTICLE_RADIUS) * (2 * PARTICLE_RADIUS);

            if (distance_sq < radius_sum_sq && distance_sq > 0) {
                float distance = std::sqrt(distance_sq);

                // ���� ��������
                float pressure_force = (particle.pressure + other->pressure) / 2 * (1 - distance / (2 * PARTICLE_RADIUS));
                dx += pressure_force * dx_diff / distance;
                dy += pressure_force * dy_diff / distance;
                dz += pressure_force * dz_diff / distance; // ���� �� z

                // ���� ��������
                float vx_diff = other->vx - particle.vx;
                float vy_diff = other->vy - particle.vy;
                float vz_diff = other->vz - particle.vz; // ������� �� vz
                float viscosity_force = VISCOSITY * (vx_diff * dx_diff + vy_diff * dy_diff + vz_diff * dz_diff) / distance_sq; // 3D ��������
                dx += viscosity_force * dx_diff;
                dy += viscosity_force * dy_diff;
                dz += viscosity_force * dz_diff; // ���� �� z

                // ���� ���������
                float cohesion_force = COHESION_STRENGTH / distance;
                dx -= cohesion_force * dx_diff;
                dy -= cohesion_force * dy_diff;
                dz -= cohesion_force * dz_diff; // ���� �� z

                // ���� ������� (��� �������� ��� 3D)
                particle.apply_spring_force(other);
            }
        }

        // ���������� ���
        particle.vx += dx / particle.density;
        particle.vy += dy / particle.density;
        particle.vz += dz / particle.density; // ���������� ���� �� z
    }
}

// ������� ��� �������������� �������� ��������� � ������� � ������ �����������
void screenToWorld(int screenX, int screenY, float depth, float& worldX, float& worldY, float& worldZ) {
    // �������� �������� �����������
    float perspectiveScale = 500.0f / (500.0f - depth);

    worldX = (screenX - WIDTH / 2) / perspectiveScale;
    worldY = (screenY - HEIGHT / 2) / perspectiveScale;
    worldZ = depth;

    // �������� ������� ������ ��� Y
    float tempX = worldX;
    worldX = tempX * cos(-cameraYaw) - worldZ * sin(-cameraYaw);
    worldZ = tempX * sin(-cameraYaw) + worldZ * cos(-cameraYaw);

    // ���������� ������� ������
    worldX += cameraX;
    worldY += cameraY;
    worldZ += cameraZ;
}

int main(int argc, char* argv[]) {
    // ������������� SDL 
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cerr << "SDL_Init Error: " << SDL_GetError() << std::endl;
        return 1;
    }

    // �������� ���� 
    SDL_Window* window = SDL_CreateWindow("3D Fluid Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        WIDTH, HEIGHT, SDL_WINDOW_SHOWN);
    if (window == nullptr) {
        std::cerr << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return 1;
    }

    // �������� ���������
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == nullptr) {
        std::cerr << "SDL_CreateRenderer Error: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }



    // �������� ���������
    SDL_Color color1 = { 0, 0, 255, 255 };
    SDL_Color color2 = { 255, 0, 0, 255 };
    create_gradient(color1, color2, GRADIENT_STEPS, color_gradient);

    // �������� ������ ������
    std::vector<Particle> particles;

    // ����� ��� ������ ����
    bool mouse_pressed_left = false;
    bool mouse_pressed_right = false;

    // ������ ����������� ������ 
    std::vector<Particle*> grabbed_particles;

    // ��������� ��������� ����� 
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distrib(-PARTICLE_RADIUS * 2, PARTICLE_RADIUS * 2);
    std::uniform_real_distribution<> z_distrib(PAD_Z, 400 - PAD_Z); // ������������� �� z

    // ���������� �������
    int cursorX, cursorY;

    // ������� ����
    bool running = true;
    Uint32 last_time = SDL_GetTicks();
    while (running) {
        // ��������� ������� 
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
                    SDL_GetMouseState(&cursorX, &cursorY);

                    // ������ ������ � ������ �����������
                    for (Particle& particle : particles) {
                        float worldX, worldY, worldZ;
                        screenToWorld(cursorX, cursorY, particle.z, worldX, worldY, worldZ);

                        float dx = particle.x - worldX;
                        float dy = particle.y - worldY;
                        float dz = particle.z - worldZ;
                        float distance_sq = dx * dx + dy * dy + dz * dz;
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
            else if (event.type == SDL_MOUSEWHEEL) {
                // �����������/��������� � ������� �������� ����
                SDL_GetMouseState(&cursorX, &cursorY);
                float worldX, worldY, worldZ;
                screenToWorld(cursorX, cursorY, 0.0f, worldX, worldY, worldZ);

                float directionX = worldX - cameraX;
                float directionY = worldY - cameraY;
                float directionZ = worldZ - cameraZ;
                float length = sqrt(directionX * directionX + directionY * directionY + directionZ * directionZ);

                if (event.wheel.y > 0) { // ��������� �����
                    cameraX += directionX / length * 10;
                    cameraY += directionY / length * 10;
                    cameraZ += directionZ / length * 10;
                }
                else if (event.wheel.y < 0) { // ��������� ����
                    cameraX -= directionX / length * 10;
                    cameraY -= directionY / length * 10;
                    cameraZ -= directionZ / length * 10;
                }
            }
            else if (event.type == SDL_KEYDOWN) {
                if (event.key.keysym.sym == SDLK_r) { // �������� ����� �� ������
                    particles.clear();
                }
            }
        }

        // ��������� ��������� ������ (��� �������� �������� ������)
        const Uint8* keystate = SDL_GetKeyboardState(NULL);
        if (keystate[SDL_SCANCODE_W]) { // ������
            cameraX += 10 * sin(cameraYaw);
            cameraZ += 10 * cos(cameraYaw);
        }
        if (keystate[SDL_SCANCODE_S]) { // �����
            cameraX -= 10 * sin(cameraYaw);
            cameraZ -= 10 * cos(cameraYaw);
        }
        if (keystate[SDL_SCANCODE_A]) { // �����
            cameraX -= 10 * cos(cameraYaw);
            cameraZ += 10 * sin(cameraYaw);
        }
        if (keystate[SDL_SCANCODE_D]) { // ������
            cameraX += 10 * cos(cameraYaw);
            cameraZ -= 10 * sin(cameraYaw);
        }
        if (keystate[SDL_SCANCODE_SPACE]) { // �����
            cameraY -= 10;
        }
        if (keystate[SDL_SCANCODE_LSHIFT]) { // ����
            cameraY += 10;
        }

        // �������� ������ � ������� Ctrl + �������� ����
        if (keystate[SDL_SCANCODE_LCTRL]) {
            SDL_SetRelativeMouseMode(SDL_TRUE); // �������� ������������� ����� ����
            int mouseX, mouseY;
            SDL_GetRelativeMouseState(&mouseX, &mouseY);
            cameraYaw += mouseX * 0.01f;
            cameraPitch += mouseY * 0.01f;
        }
        else {
            SDL_SetRelativeMouseMode(SDL_FALSE); // ��������� ������������� ����� ����
        }

        // �������� ����� ������ (��������� z ����������)
        if (mouse_pressed_left && particles.size() < MAX_PARTICLES) {
            for (int i = 0; i < PARTICLE_CREATION_RATE; ++i) {
                int mouseX, mouseY;
                SDL_GetMouseState(&mouseX, &mouseY);
                float worldX, worldY, worldZ;
                screenToWorld(mouseX, mouseY, 200.0f, worldX, worldY, worldZ); // ���������� 200.0f ��� worldZ
                particles.emplace_back(worldX + distrib(gen), worldY + distrib(gen), worldZ);
            }
        }

        // ����������� ����������� ������ � ������ �����������
        if (mouse_pressed_right) {
            SDL_GetMouseState(&cursorX, &cursorY);
            for (Particle* particle : grabbed_particles) {
                float worldX, worldY, worldZ;
                screenToWorld(cursorX, cursorY, particle->z, worldX, worldY, worldZ);

                float dx = worldX - particle->x;
                float dy = worldY - particle->y;
                float dz = worldZ - particle->z;
                particle->vx = dx * DRAG_COEFFICIENT;
                particle->vy = dy * DRAG_COEFFICIENT;
                particle->vz = dz * DRAG_COEFFICIENT;
            }
        }

        // ���������� ������
        find_neighbors(particles);
        calculate_density_pressure(particles);
        calculate_forces(particles);

        // ���������� ������
        for (Particle& particle : particles) {
            particle.update();
        }

        // ���������
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // ������ ����������� (�����)
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        // ��������� ����� (�� �����������)
        // ...

        // ������ ������� 
        for (Particle& particle : particles) {
            particle.draw(renderer);
        }

        // ���������� �������
        int cursorX = 0, cursorY = 0; // �������������� ��� ����������
        // ������ ����� ����� ��� �������
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderDrawPoint(renderer, cursorX, cursorY);



        // ����� ���������� � ��������� ����
        SDL_GetMouseState(&cursorX, &cursorY);
        float worldX, worldY, worldZ;
        screenToWorld(cursorX, cursorY, 0.0f, worldX, worldY, worldZ);

        std::string title = "3D Fluid Simulation | Particles: " + std::to_string(particles.size()) + "/" + std::to_string(MAX_PARTICLES) +
            " | Camera: (" + std::to_string(int(cameraX)) + ", " + std::to_string(int(cameraY)) + ", " + std::to_string(int(cameraZ)) +
            ") Yaw: " + std::to_string(int(cameraYaw * 180.0f / M_PI)) + "� Pitch: " + std::to_string(int(cameraPitch * 180.0f / M_PI)) + "�" +
            " | Cursor: (" + std::to_string(int(worldX)) + ", " + std::to_string(int(worldY)) + ", " + std::to_string(int(worldZ)) + ")";

        SDL_SetWindowTitle(window, title.c_str());

        // ����� FPS � ���������� ������ 
        Uint32 current_time = SDL_GetTicks();
        float fps = 1000.0f / (current_time - last_time);
        last_time = current_time;

        SDL_RenderPresent(renderer);
        SDL_Delay(1000 / FPS);
    }

    // ������������ �������� 
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}