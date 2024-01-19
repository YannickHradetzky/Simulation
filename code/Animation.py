import os
import numpy as np
import matplotlib.pyplot as plt
import imageio

def read_data(positions_dir, velocities_dir):
    position_files = sorted(os.listdir(positions_dir), key=lambda x: int(x.split('.')[0]))
    velocity_files = sorted(os.listdir(velocities_dir), key=lambda x: int(x.split('.')[0]))

    positions = [np.loadtxt(os.path.join(positions_dir, f)) for f in position_files]
    velocities = [np.loadtxt(os.path.join(velocities_dir, f)) for f in velocity_files]

    return positions, velocities

def find_limits(positions):
    x_min = np.min([np.min(pos[:, 0]) for pos in positions])
    x_max = np.max([np.max(pos[:, 0]) for pos in positions])
    y_min = np.min([np.min(pos[:, 1]) for pos in positions])
    y_max = np.max([np.max(pos[:, 1]) for pos in positions])

    return x_min, x_max, y_min, y_max

def create_particle_animation(positions_dir, velocities_dir, output_gif):
    positions, velocities = read_data(positions_dir, velocities_dir)
    x_min, x_max, y_min, y_max = find_limits(positions)

    # Create a figure
    fig, ax = plt.subplots()
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    with imageio.get_writer(output_gif, duration=0.1) as writer:
        for frame in range(len(positions)):
            # Print progress
            print(f'Creating frame {frame}')

            ax.cla()
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)
            ax.set_title(f'Frame {frame}')

            # Load position and velocity data for each particle
            pos_data = positions[frame]
            vel_data = velocities[frame]

            # Extract the position and velocity for the current frame
            x = pos_data[:, 0]
            y = pos_data[:, 1]
            u = vel_data[:, 0]
            v = vel_data[:, 1]
            ax.quiver(x, y, u, v)

            # Save each frame to the GIF
            frame_filename = f"frame_{frame:04d}.png"
            plt.savefig(frame_filename)  # Save the figure
            writer.append_data(imageio.imread(frame_filename))
            os.remove(frame_filename)  # Remove the temporary file


    # Clean up temporary image files
    for frame in range(len(positions)):
        os.remove(f"frame_{frame:04d}.png")

if __name__ == "__main__":
    pos_dir = '/Users/yannick/Documents/Master/Many Particle Physics/Simulation/out/vicsek/Positions/N1000/2.000000'
    vel_dir = '/Users/yannick/Documents/Master/Many Particle Physics/Simulation/out/vicsek/Velocities/N1000/2.000000'
    output_gif_file = '/Users/yannick/Documents/Master/Many Particle Physics/Simulation/out/vicsek/2.gif'

    create_particle_animation(pos_dir, vel_dir, output_gif_file)
