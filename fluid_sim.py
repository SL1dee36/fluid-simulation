# VERSION 0.7
# Created by @Sl1dee36
# Simple Fluid Simulation

import pygame
import random

# Настройки
WIDTH = 800
HEIGHT = 600
FPS = 60
PARTICLE_RADIUS = 4
PARTICLE_COLOR = (0, 255, 255)
GRAVITY = 0.3
REST_DENSITY = 20  # Целевая плотность жидкости
GAS_CONSTANT = 1  # Сила, отталкивающая частицы друг от друга
VISCOSITY = 25  # Сопротивление движению жидкости
PARTICLE_CREATION_RATE = 15  # Количество частиц, создаваемых за кадр
COHESION_STRENGTH = 0.001  # Сила сцепления
DAMPING = 0.99  # Коэффициент затухания скорости
DRAG_COEFFICIENT = 0.1  # Коэффициент сопротивления при захвате
GRAB_RADIUS = 50  # Радиус захвата в пикселях
SPRING_CONSTANT = 0.05  # Жесткость пружины для soft-body

# Препятствие
OBSTACLE_X = 400
OBSTACLE_Y = 300
OBSTACLE_RADIUS = 100

# Стены
WALL_THICKNESS = 0
WALL_COLOR = (255, 255, 255)
PAD_X = 20
PAD_Y = 20

# Инициализация Pygame
pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("2D Fluid Simulation")
clock = pygame.time.Clock()

# Функция для создания градиента цветов
def create_gradient(color1, color2, steps):
    gradient = []
    for i in range(steps):
        r = int(color1[0] * (1 - i / (steps - 1)) + color2[0] * i / (steps - 1))
        g = int(color1[1] * (1 - i / (steps - 1)) + color2[1] * i / (steps - 1))
        b = int(color1[2] * (1 - i / (steps - 1)) + color2[2] * i / (steps - 1))
        gradient.append((r, g, b))
    return gradient

# Градиент цветов от синего к красному
color_gradient = create_gradient((0, 0, 255), (255, 0, 0), 256)
max_speed = 20  # Максимальная скорость для расчета цвета


# Класс частицы
class Particle:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.vx = 0
        self.vy = 0
        self.density = 0
        self.pressure = 0
        self.near_particles = []
        self.grabbed = False

    def update(self):
        self.vy += GRAVITY
        self.x += self.vx
        self.y += self.vy

        # Проверка столкновений с границами экрана
        if self.x < PARTICLE_RADIUS + PAD_X:
            self.x = PARTICLE_RADIUS + PAD_X
            self.vx *= -0.7
        elif self.x > WIDTH - PARTICLE_RADIUS - PAD_X:
            self.x = WIDTH - PARTICLE_RADIUS - PAD_X
            self.vx *= -0.7

        if self.y > HEIGHT - PARTICLE_RADIUS - PAD_Y:
            self.y = HEIGHT - PARTICLE_RADIUS - PAD_Y
            self.vy *= -0.7

        # Проверка столкновения с препятствием
        distance_to_obstacle = ((self.x - OBSTACLE_X) ** 2 + (self.y - OBSTACLE_Y) ** 2) ** 0.5
        if distance_to_obstacle < PARTICLE_RADIUS + OBSTACLE_RADIUS:
            overlap = (PARTICLE_RADIUS + OBSTACLE_RADIUS) - distance_to_obstacle
            self.x += overlap * (self.x - OBSTACLE_X) / distance_to_obstacle
            self.y += overlap * (self.y - OBSTACLE_Y) / distance_to_obstacle
            self.vx *= -0.7
            self.vy *= -0.7

        # Проверка на пересечение с другими частицами
        for other in particles:
            if self != other:
                distance = ((self.x - other.x) ** 2 + (self.y - other.y) ** 2) ** 0.5
                if distance < 0.1:
                    distance == 0.1
                if distance < 2 * PARTICLE_RADIUS and distance != 0:
                    # Частицы пересекаются, "разъединяем" их
                    overlap = 2 * PARTICLE_RADIUS - distance
                    self.x += overlap / 2 * (self.x - other.x) / distance
                    self.y += overlap / 2 * (self.y - other.y) / distance
                    other.x -= overlap / 2 * (self.x - other.x) / distance
                    other.y -= overlap / 2 * (self.y - other.y) / distance

        # Демпфирование скорости
        self.vx *= DAMPING
        self.vy *= DAMPING

    def draw(self):
        # Использование градиента для плавного перехода цветов
        speed = (self.vx ** 2 + self.vy ** 2) ** 0.5
        color_index = min(int(speed / max_speed * (len(color_gradient) - 1)), len(color_gradient) - 1)
        color = color_gradient[color_index]
        pygame.draw.circle(screen, color, (int(self.x), int(self.y)), PARTICLE_RADIUS)

    def reset_velocity(self):
        self.vx = 0
        self.vy = 0

    def apply_spring_force(self, other):
        distance = ((self.x - other.x) ** 2 + (self.y - other.y) ** 2) ** 0.5
        if distance > 2 * PARTICLE_RADIUS:
            force = SPRING_CONSTANT * (distance - 2 * PARTICLE_RADIUS)
            dx = force * (self.x - other.x) / distance
            dy = force * (self.y - other.y) / distance
            self.vx -= dx / self.density
            self.vy -= dy / self.density
            other.vx += dx / other.density
            other.vy += dy / other.density

# Функция для нахождения ближайших частиц (оптимизация)
def find_neighbors(particles):
    for i in range(len(particles)):
        particles[i].near_particles = []
        for j in range(i + 1, len(particles)):
            distance = ((particles[i].x - particles[j].x) ** 2 + (particles[i].y - particles[j].y) ** 2) ** 0.5
            if distance < 2 * PARTICLE_RADIUS:
                particles[i].near_particles.append(particles[j])
                particles[j].near_particles.append(particles[i])


# Функция для вычисления плотности и давления
def calculate_density_pressure(particles):
    for i in range(len(particles)):
        particles[i].density = 0
        for j in range(len(particles[i].near_particles)):
            distance = ((particles[i].x - particles[j].x) ** 2 + (particles[i].y - particles[j].y) ** 2) ** 0.5
            if distance < 2 * PARTICLE_RADIUS:
                particles[i].density += 1  # Упрощенное вычисление плотности

        # Добавляем минимальную плотность, чтобы избежать деления на ноль
        particles[i].density = max(particles[i].density, 0.1)
        particles[i].pressure = GAS_CONSTANT * (particles[i].density - REST_DENSITY)


# Функция для вычисления сил взаимодействия
def calculate_forces(particles):
    for i in range(len(particles)):
        dx = dy = 0
        for j in range(len(particles[i].near_particles)):
            distance = ((particles[i].x - particles[j].x) ** 2 + (particles[i].y - particles[j].y) ** 2) ** 0.5
            if distance < 2 * PARTICLE_RADIUS and distance > 0:
                # Вычисление силы давления
                pressure_force = (particles[i].pressure + particles[j].pressure) / 2 * (1 - distance / (2 * PARTICLE_RADIUS))
                dx += pressure_force * (particles[i].x - particles[j].x) / distance
                dy += pressure_force * (particles[i].y - particles[j].y) / distance

                # Вычисление силы вязкости
                vx_diff = particles[j].vx - particles[i].vx
                vy_diff = particles[j].vy - particles[i].vy
                viscosity_force = VISCOSITY * (vx_diff * (particles[i].x - particles[j].x) + vy_diff * (
                            particles[i].y - particles[j].y)) / distance ** 2
                dx += viscosity_force * (particles[i].x - particles[j].x)
                dy += viscosity_force * (particles[i].y - particles[j].y)

                # Вычисление силы сцепления (cohesion)
                cohesion_force = COHESION_STRENGTH / distance
                dx -= cohesion_force * (particles[i].x - particles[j].x)
                dy -= cohesion_force * (particles[i].y - particles[j].y)

                # Применение силы пружины (soft-body)
                particles[i].apply_spring_force(particles[j])

        # Применение сил к скорости частицы
        particles[i].vx += dx / particles[i].density
        particles[i].vy += dy / particles[i].density

# Создание списка частиц
particles = []

# Главный цикл
running = True
mouse_pressed_left = False
mouse_pressed_right = False
grabbed_particles = []
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        elif event.type == pygame.MOUSEBUTTONDOWN:
            if event.button == 1:  # Левая кнопка мыши
                mouse_pressed_left = True
            elif event.button == 3:  # Правая кнопка мыши
                mouse_pressed_right = True
                mouse_x, mouse_y = pygame.mouse.get_pos()
                for particle in particles:
                    distance_to_mouse = ((particle.x - mouse_x)**2 + (particle.y - mouse_y)**2)**0.5
                    if distance_to_mouse <= GRAB_RADIUS:
                        particle.grabbed = True
                        grabbed_particles.append(particle)
        elif event.type == pygame.MOUSEBUTTONUP:
            if event.button == 1:
                mouse_pressed_left = False
            elif event.button == 3:
                mouse_pressed_right = False
                for particle in grabbed_particles:
                    particle.grabbed = False
                grabbed_particles = []
        elif event.type == pygame.KEYDOWN:
            if event.key == pygame.K_s:  # Нажата клавиша "S"
                for particle in particles:
                    particle.reset_velocity() 

    # Создание новых частиц при нажатой левой кнопке мыши
    if mouse_pressed_left:
        for _ in range(PARTICLE_CREATION_RATE):
            mouse_x, mouse_y = pygame.mouse.get_pos()
            x = random.randint(mouse_x - PARTICLE_RADIUS * 2, mouse_x + PARTICLE_RADIUS * 2)
            y = random.randint(mouse_y - PARTICLE_RADIUS * 2, mouse_y + PARTICLE_RADIUS * 2)
            particles.append(Particle(x, y))

    # Перемещение захваченных частиц
    if mouse_pressed_right:
        mouse_x, mouse_y = pygame.mouse.get_pos()
        for particle in grabbed_particles:
            dx = mouse_x - particle.x
            dy = mouse_y - particle.y
            particle.vx = dx * DRAG_COEFFICIENT
            particle.vy = dy * DRAG_COEFFICIENT

    # Обновление физики жидкости
    find_neighbors(particles)
    calculate_density_pressure(particles)
    calculate_forces(particles)

    # Обновление частиц
    for particle in particles:
        particle.update()

    # Отрисовка
    screen.fill((0, 0, 0))

    # Рисуем стены
    pygame.draw.rect(screen, WALL_COLOR, (PAD_X, PAD_Y, WIDTH - 2 * PAD_X, WALL_THICKNESS))  # Верх
    pygame.draw.rect(screen, WALL_COLOR, (PAD_X, HEIGHT - PAD_Y - WALL_THICKNESS, WIDTH - 2 * PAD_X, WALL_THICKNESS))  # Низ
    pygame.draw.rect(screen, WALL_COLOR, (PAD_X, PAD_Y, WALL_THICKNESS, HEIGHT - 2 * PAD_Y))  # Лево
    pygame.draw.rect(screen, WALL_COLOR, (WIDTH - PAD_X - WALL_THICKNESS, PAD_Y, WALL_THICKNESS, HEIGHT - 2 * PAD_Y))  # Право

    # Рисуем препятствие
    pygame.draw.circle(screen, (255, 255, 255), (OBSTACLE_X, OBSTACLE_Y), OBSTACLE_RADIUS)

    # Рисуем частицы 
    for particle in particles:
        particle.draw()
    
    # Вывод FPS в заголовке окна
    fps = clock.get_fps()
    pygame.display.set_caption(f"2D Fluid Simulation - FPS: {fps:.2f}")
    
    pygame.display.flip()
    clock.tick(FPS)

pygame.quit()