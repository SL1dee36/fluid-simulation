# 2D Fluid Simulation

**Version:** 0.975 Raw

**Author:** [SL1dee36](https://github.com/SL1dee36)

This program simulates a 2D fluid using the Smoothed Particle Hydrodynamics (SPH) method. It allows users to interact with the fluid by creating new particles and manipulating existing ones.

Currect verison demo:
https://www.youtube.com/watch?v=xbR4eWPZM4s

Old verison demo:
https://github.com/user-attachments/assets/170b9447-a075-4b4f-9c04-e5bc13744b59

## Features

* **Realistic fluid behavior:** Simulates fluid movement, pressure, viscosity, and interaction with boundaries.
* **Interactive controls:**
    * **Left Mouse Button:** Creates new particles at the cursor position.
    * **Right Mouse Button:** Drags and interacts with particles within a radius.
    * **'S' Key:** Resets the velocity of all particles.
* **Obstacle interaction:** Includes a circular obstacle that particles can collide with.
* **Visual representation:** Renders particles with a color gradient based on their speed.
* **Performance optimizations:** Implements efficient neighbor search and force calculation algorithms.

## Dependencies

* **SDL2:** Simple DirectMedia Layer 2 library for graphics rendering and input handling.

## Building and Running

1. **Install SDL2:** Instructions for installing SDL2 can be found on the official website: [https://www.libsdl.org/download-2.0.php](https://www.libsdl.org/download-2.0.php)
2. **Compile the code:** Use a C++ compiler that supports SDL2 (e.g., g++):
   ```bash
   g++ -o fluid_simulation fluid_simulation.cpp -lSDL2
   ```
3. **Run the executable:**
   ```bash
   ./fluid_simulation
   ```

## Configuration

The simulation parameters can be adjusted by modifying the constants at the beginning of the code:

* **`WIDTH` and `HEIGHT`:** Dimensions of the simulation window.
* **`FPS`:** Frames per second.
* **`PARTICLE_RADIUS`:** Radius of each particle.
* **`GRAVITY`:** Strength of the gravitational force.
* **`REST_DENSITY`:** Rest density of the fluid.
* **`GAS_CONSTANT`:** Gas constant used in pressure calculation.
* **`VISCOSITY`:** Viscosity of the fluid.
* **`PARTICLE_CREATION_RATE`:** Number of particles created per frame when the left mouse button is held down.
* **`COHESION_STRENGTH`:** Strength of the cohesion force between particles.
* **`DAMPING`:** Damping factor applied to particle velocities.
* **`DRAG_COEFFICIENT`:** Drag coefficient applied to grabbed particles.
* **`GRAB_RADIUS`:** Radius within which particles are affected by the right mouse button.
* **`SPRING_CONSTANT`:** Spring constant used in the spring force calculation.
* **`OBSTACLE_X` and `OBSTACLE_Y`:** Position of the circular obstacle.
* **`OBSTACLE_RADIUS`:** Radius of the circular obstacle.
* **`WALL_THICKNESS`:** Thickness of the walls (currently set to 0 for optimized rendering).
* **`WALL_COLOR`:** Color of the walls.
* **`PAD_X` and `PAD_Y`:** Padding between the walls and the simulation area.
* **`GRADIENT_STEPS`:** Number of steps in the color gradient.
* **`max_speed`:** Maximum speed for color mapping.

## Known Issues

* Performance may degrade with a large number of particles.
* The simulation is not perfectly accurate and may exhibit some unrealistic behavior in certain scenarios.

## Future Improvements

* Implement more advanced SPH features, such as surface tension and boundary handling.
* Optimize performance further to handle more particles.
* Add more interactive features and customization options.

## License

This project is released under the MIT License. See the LICENSE file for details.
