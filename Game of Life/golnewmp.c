#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <omp.h>

void initialWorld(int* world_, const int height, const int width);
int countLivingNeighbours(int*  world_, const int width, int cell_row, int cell_col);
void computeNextGeneration(int* world_, const int height, const int width);
void copyWorld(int* __restrict old, int* __restrict new, const int height, const int width);
void printWorld(int* world_, const int height, const int width);

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
} 

int main(int argc, char *argv[]){
    double time1 = get_wall_seconds();
    if (argc != 4){
        printf("Check your input..\n");
        return -1;
    }
    const int N = atoi(argv[1]);
    const int HEIGHT = atoi(argv[2])+2;
    const int WIDTH = atoi(argv[3])+2;
    int *world = (int*)calloc(HEIGHT*WIDTH,sizeof(int*));
    initialWorld(world, HEIGHT, WIDTH);
    //printf("\n Intial generation:\n");
    //printWorld(world, HEIGHT, WIDTH); 
    for(int count = 1; count<N; count++){
        //printf("\n Generation number: %d\n", count);
        computeNextGeneration(world, HEIGHT, WIDTH);
        //printWorld(world, HEIGHT, WIDTH);
    }
    printf(" The game of life main took %7.3f wall seconds.\n", get_wall_seconds()-time1);
    free(world);
    
    return 0;
}

/**
 * Randomly generates the initial board with 0's and 1's,
 * with size rows*columns
 * 0 = dead, 1 = alive
 */ 
void initialWorld(int* world_, const int height, const int width){
    #pragma omp parallel for
	for (int row = 1; row < height-1; row++){
		for (int col = 1; col < width-1; col++){
         	world_[col+width*row] = rand() % 2; 
        }
    }     
}

/**
 * Counting the number of living neigbour of a given cell,
 * with position world[cell_row][cell_col]
 * Iterating through the 8 neighbours 
 * If the cell looking at is outside the world or its the current cell, skip it and continue
 */
int countLivingNeighbours(int* world_, const int width, int cell_row, int cell_col) {
	int neighbours = 0; 
	for (int i = -1; i < 2; i+=2) {
		for (int j = -1; j < 2; j++) {

            neighbours += world_[cell_col+j+width*(cell_row+i)];
        }
    }
    neighbours += world_[cell_col+1+width*(cell_row)];
    neighbours += world_[cell_col-1+width*(cell_row)];

	return neighbours;
}
/**
* The game of life rules: 
* If an alive cell has 2 or 3 living neighbours, it surives
* If a dead cell has 3 neighbours, it becomes a live cell
* All other cells die in the next generation
* Using function @countLivingNeighbours to count the number of living neighbours for that cell.
*/
void computeNextGeneration(int* world_, const int height, const int width){
    int *nextworld = (int*)calloc(height*width,sizeof(int*));
    #pragma omp parallel for 
    for (int row = 1; row < height-1; row++) {
		for (int col = 1; col < width-1; col++) {
			int nr_neighbours = countLivingNeighbours(world_, width, row, col);
			if (nr_neighbours == 3) {
				nextworld[col+width*row] = 1;
			} else if (world_[col+width*row] == 1 && nr_neighbours == 2) {
				nextworld[col+width*row] = 1;
			} else {
				nextworld[col+width*row] = 0;
			}
		}
	}
    copyWorld(world_, nextworld, height, width);

    free(nextworld);
}

/**
 * Copies the world, so that the new generation becomes the old generation for the next round.
*/
void copyWorld(int* __restrict old, int* __restrict new, const int height, const int width){
    #pragma omp parallel for
    for (int row = 1; row < height-1; row++) {
		for (int col = 1; col < width-1; col++) {
            old[col+width*row] = new[col+width*row];
        }
    }        
}

/**
 * Prints the world, where # = Alive, and 0 = Dead
 * 
*/
void printWorld(int* world_, const int height, const int width) {
    printf("\n\n\n");
   	for (int row = 1; row < height-1; row++) {
		for (int col = 1; col < width-1; col++) {
            if(world_[col+width*row] == 1){
			printf("#");
		}
            else{
                printf("0");
            }
        }
	printf("\n");
    }
}
