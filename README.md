<a name="readme-top"></a>

<h2 align="center">LINMA1170 : Project</h2>

</div>

## Description of the project

This program is designed to optimize the dimensions of a tuning fork based on given target parameters. The program takes command-line arguments for the results and animation filenames, if animation is desired. Many parameters can be tweaked in the main function of project.c as targetted mass, tolerance on the frequency error, ...

## Prerequisites

The project needs OpenMP for multithreading. It is assumed to be already installed on Linux machines. 

To save the animations, you will need FFmpeg:
* FFmpeg
  ```sh
  sudo apt install ffmpeg
  ```

## Usage
1. To compile the project, use the following command:

   Note: The Makefile is made for Linux systems.

  ```sh
  make all
  ```
2. To run the program, use the following command:
  ```sh
  ./project <results_filename> <animation_filename>
  ```
- results_filename: The output file to write the dimensions of the tuning fork.
- animation_filename (optional): The output file to save the animation (you need to write the format extension, for example .mp4).

Note: If the animation filename is not provided, the program will only output the tuning fork dimensions and skip the animation generation.

3. Once ran, the program saves the dimensions of the tuning fork in /results/results_filename. The animation, if generated, is saved in /animation/animation_filename.

## Note for Windows users
If the animation generation does not work for Windows users, an animation is already provided in the "animation" folder.
