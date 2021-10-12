//
//  main.c
//  mknapsack
//
//  Created by bai on 29/03/2019.
//  Copyright © 2019 UNNC. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

/* parameters */
int RAND_SEED[] = {1,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
int NUM_OF_RUNS = 1;
static int POP_SIZE = 100; //global parameters
int MAX_NUM_OF_GEN = 1000; //max number of generations
int MAX_TIME = 60;  //max amount of time permited (in sec)
float CROSSOVER_RATE = 0.8;
float MUTATION_RATE = 0.1;
struct solution_struct* best_sln;  //global best solution

struct item_struct{
    int dim; //no. of dimensions
    int *size; //volume of item in all dimensions
    int p; //profit 
    float *ratio;
};

struct problem_struct{
    int n; //number of items
    int dim; //number of dimensions
    struct item_struct* items;
    int* capacities;  //knapsack capacities
    int best_sol;
};

struct solution_struct{
    struct problem_struct* prob; //maintain a shallow copy of the problem data
    int objective;
    int feasibility; //indicate the feasiblity of the solution
    int* x; //chromosome vector
    int* cap_left; //capacity left in all dimensions
};

//return a random number between 0 and 1
float rand_01()
{
    float number;
    number = (float) rand();
    number = number/RAND_MAX;
    //printf("rand01=%f\n", number);
    return number;
}

//return a random nunber ranging from min to max (inclusive)
int rand_int(int min, int max)
{
    int div = max-min+1;
    int val =rand() % div + min;
    //printf("rand_range= %d \n", val);
    return val;
}

//a random number between 0-1
float rand_float(){
    return rand()/(RAND_MAX+1.0); 
}

void free_problem(struct problem_struct* prob)
{
    if(prob!=NULL)
    {
        if(prob->capacities !=NULL){
            free(prob->capacities);
        }
        if(prob->items!=NULL)
        {
            for(int j=0; j<prob->n; j++)
            {
                if(prob->items[j].size != NULL){
                    free(prob->items[j].size);
                }
                if(prob->items[j].ratio != NULL){
                    free(prob->items[j].ratio);
                }
            }
            free(prob->items);
        }
        free(prob);
    }
}

void free_population(struct solution_struct* pop, int size)
{
    for(int p=0; p<size; p++)
    {
        free(pop[p].x);
        free(pop[p].cap_left);
        free(pop[p].prob);
    }
}

void init_problem(int n, int dim, int best, struct problem_struct** my_prob)
{
    struct problem_struct* new_prob = malloc(sizeof(struct problem_struct));
    if (new_prob == NULL){
        printf("error 1\n");
    }
    
    new_prob->n=n; 
    new_prob->dim=dim;
    new_prob->best_sol=best;
    new_prob->items=malloc(sizeof(struct item_struct)*n);
    if (new_prob->items == NULL){
        perror("error 2");
    }
    
    for(int j=0; j<n; j++){
        new_prob->items[j].size = malloc(sizeof(int)*dim);
        new_prob->items[j].ratio = malloc(sizeof(float)*dim);
        if (new_prob->items[j].size == NULL || new_prob->items[j].ratio == NULL){
            perror("error 3");
        }
        
    }
    new_prob->capacities = malloc(sizeof(int)*dim);
    if (new_prob->capacities == NULL){
        perror("error 4");
    }
    
    *my_prob = new_prob;
}

struct problem_struct** load_problems(int num_of_probs, FILE *filepointer)
{
    int i,j,k, n, dim, best, profit, volume, cap;
    struct problem_struct** my_problems = malloc(sizeof(struct problem_struct)*num_of_probs);
    if (my_problems == NULL){ 
        printf("exit 1\n");
        exit(1);
    }
    for(k=0; k<num_of_probs; k++)
    {
        fscanf(filepointer, "%d", &n);
        fscanf(filepointer, "%d", &dim);
        fscanf(filepointer, "%d", &best);
        init_problem(n, dim, best, &my_problems[k]); 
        //printf("number of items: %d, dim: %d, best: %d\n", my_problems[k]->n, my_problems[k]->dim, my_problems[k]->best_sol);
        
        //for every item in problem k, initialize profit,
        for(j=0; j<n; j++)
        {
            fscanf(filepointer, "%d", &profit);
            my_problems[k]->items[j].dim=dim;
            my_problems[k]->items[j].p=profit;
            //printf("item %d, dim %d, profit %d\n", j, my_problems[k]->items[j].dim, my_problems[k]->items[j].p);
        }
        //for every item in problem k, initialize volume and ratio in each dimension 
        for (i=0; i < dim; i++)
        {
            for(j=0; j<n; j++)
            {
                fscanf(filepointer, "%d", &volume);
                my_problems[k]->items[j].size[i] = volume;
                my_problems[k]->items[j].ratio[i] = (float)my_problems[k]->items[j].p/(float)volume;
                //printf("volume %d, ratio %f\n", my_problems[k]->items[j].size[i], my_problems[k]->items[j].ratio[i]);
            }
        }
        //for every dimension, initialize capacity 
        for(i=0; i < dim; i++){
            fscanf(filepointer, "%d", &cap);
            my_problems[k]->capacities[i] = cap;
            //printf("prob %d, dim %d, capacity: %d\n", k, my_problems[k]->dim, my_problems[k]->capacities[i]);
        }
    }
    fclose(filepointer);
    return my_problems;
}

void evaluate_solution(struct solution_struct* sln)
{
    //evaluate the feasibility and objective of the solution
    sln->objective = 0; 
    sln->feasibility = 1;
    struct item_struct* items_p = sln->prob->items;
    
    for(int i=0; i< items_p->dim; i++)
    {
        sln->cap_left[i] = sln->prob->capacities[i];
        for(int j=0; j < sln->prob->n; j++)
        {
            sln->cap_left[i] -= items_p[j].size[i]*sln->x[j];
            if(sln->cap_left[i]<0) {
                sln->feasibility = -1*i; //exceeding capacity
                return;
            }
        }
    }
    if(sln->feasibility > 0)
    {
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->objective += sln->x[j] * items_p[j].p;
        }
    }
}

//output a given solution to a file
void output_solution(struct solution_struct* sln, char* filename)
{
    FILE* filepointer = fopen(filename, "a");
    fprintf(filepointer, "%d\n", sln->objective);
    for (int i = 0; i < sln->prob->n; i++)
    {
        fprintf(filepointer, "%d ", sln->x[i]);
    }
    fprintf(filepointer, "\n");
    fclose(filepointer);
}

//intialise the population with random solutions
void init_population(struct problem_struct* prob, struct solution_struct* pop)
{
    for(int p=0; p<POP_SIZE; p++)
    {
        pop[p].prob = prob;
        pop[p].x = malloc(sizeof(int)*prob->n);
        pop[p].cap_left = malloc(sizeof(int)*prob->dim);
        for(int j=0; j<prob->n; j++)    
            pop[p].x[j] = 0;
        for(int i=0; i<prob->dim; i++)  
            pop[p].cap_left[i]=prob->capacities[i];
        /* create a random initial x that is feasible */
        int j=rand_int(0, prob->n-1);
        while(true)
        {
            while(pop[p].x[j]==1)
            {
                j=rand_int(0, prob->n-1); //select an unpacked item randomly
            }
            //printf("trying item %d to pack. \n", j);
            pop[p].x[j]=1;
            bool can_pack=true;
            for(int i=0; i< prob->dim; i++)
            {
                pop[p].cap_left[i] -= prob->items[j].size[i];
                if(pop[p].cap_left[i] <0) 
                    can_pack=false;
            }
            if(!can_pack)
            {   //unpack item i
                //printf("packing item %d failed. random initialisation stoped.\n", j);
                pop[p].x[j]=0;
                for(int i=0; i< prob->dim; i++)
                    pop[p].cap_left[i] += prob->items[j].size[i];
                break;
            }
        }
        evaluate_solution (&pop[p]);
    }
    
    
    for(int p=0; p<POP_SIZE; p++)
    {
        output_solution(&pop[p], "solution_file.txt");
    }
    
}
//copy a solution from another solution
bool copy_solution(struct solution_struct* dest_sln, struct solution_struct* source_sln)
{
    if(source_sln ==NULL) return false;
    if(dest_sln==NULL)
    {
        dest_sln = malloc(sizeof(struct solution_struct));
    }
    else{
        free(dest_sln->cap_left);
        free(dest_sln->x);
    }
    int n = source_sln->prob->n;
    int m =source_sln->prob->dim;
    dest_sln->x = malloc(sizeof(int)*n);
    dest_sln->cap_left=malloc(sizeof(int)*m);
    for(int i=0; i<m; i++)
        dest_sln->cap_left[i]= source_sln->cap_left[i];
    for(int j=0; j<n; j++)
        dest_sln->x[j] = source_sln->x[j];
    dest_sln->prob= source_sln->prob;
    dest_sln->feasibility=source_sln->feasibility;
    dest_sln->objective=source_sln->objective;
    return true;
}

//select individuals with higher fitness 
void selection(struct solution_struct* curt_pop)
{
    //calculate fitness for each individual
    int *roulette= malloc(sizeof(int)*POP_SIZE);
    //temperorily store the selected population 
    struct solution_struct* new_pop = malloc(sizeof(struct solution_struct)*POP_SIZE);

    for (int p = 0; p < POP_SIZE; p++){
        copy_solution(&new_pop[p], &curt_pop[p]);
    }
    roulette[0] = curt_pop[0].objective;

    for(int p = 1; p < POP_SIZE; p++){
        roulette[p] = roulette[p-1] + curt_pop[p].objective;
        //printf("roulette %d: %d\n",p-1, roulette[p-1]);
    }
    int totalFitness = roulette[POP_SIZE-1];

    //apply roulette wheel selection  
    for (int p = 0; p < POP_SIZE; p++)
    {
        int dart = rand_int(0, totalFitness-1), q = 0;
        while (dart >= roulette[q] && q < POP_SIZE)
        {
            if (dart < roulette[q+1]){
                //printf("-------->> dart = %d\n", dart);
                //printf("cur pop = %d\n", q);
                copy_solution(&new_pop[p], &curt_pop[q]);
                break;
            }
            q++;
        }
    }

    for (int p = 0; p < POP_SIZE; p++){
        copy_solution(&curt_pop[p], &new_pop[p]);
        //printf("%d, %d, \n", curt_pop[p].objective, new_pop[p].objective);
    }
    free(roulette);
    free(new_pop);
}
//generate a new population
void cross_over(struct solution_struct* new_pop, struct solution_struct* curt_pop)
{
    //pairing and cross-over
    int pool = POP_SIZE-1, first, second, xpt;
    int *temp = malloc(sizeof(int)*curt_pop->prob->n); //temperorily store a piece of chromosome;

    while (pool >= 0){
        first = rand_int(0, pool);    //randomly pick 2 for mating
        second = rand_int(0, pool);
        
        float probability = rand_float();
        if(probability > CROSSOVER_RATE){
            pool -= 2;
            continue;
        }

        //generate random cross-over point and exchange 
        xpt = rand_int(0, curt_pop->prob->n-1);
        for (int i = 0; i < xpt; i++)
        {
            temp[i] = curt_pop[first].x[i];
            curt_pop[first].x[i] = curt_pop[second].x[i];
            curt_pop[second].x[i] = temp[i];
        }
        pool -= 2;
    }
    for (int i = 0; i < POP_SIZE; i++){
        copy_solution(&new_pop[i], &curt_pop[i]);
    }
    free(temp);
}

//apply mutation to a population
void mutation(struct solution_struct* pop)
{
    MUTATION_RATE = (float)1/(float)pop->prob->n;
    //printf("mutation rate: %f\n", MUTATION_RATE);
    float *probability = malloc(sizeof(float)*POP_SIZE);
    int mutation_point;

    for (int p = 0; p < POP_SIZE; p++){
        probability[p] = rand_float();
        //printf("pop: %d, mutation rate: %f\n", p, probability[p]);
        //if mutate, randomly pick an allele for mutation
        if (probability[p] < MUTATION_RATE){
            mutation_point = rand_int(0, pop->prob->n-1);
            //printf("%d\n", mutation_point);
            if (pop[p].x[mutation_point] == 1)
                pop[p].x[mutation_point] = 0;
            else
                pop[p].x[mutation_point] = 1;
        }
    }
    free(probability);
}
//find the items in the knapsack with smallest ratio in all dimensions 
int min_ratio(struct solution_struct* solution)
{
    float cur = 0, min = 0;
    int index;
    for (int d = 0; d < solution->prob->dim; d++)
        min += solution->prob->items[0].ratio[d];  //initialize min
    
    for (int i = 1; i < solution->prob->n; i++){
        for (int d = 0; d < solution->prob->dim; d++){
            cur += solution->prob->items[i].ratio[d];
        } 
        if (min > cur && solution->x[i]== 1){
            min = cur;
            index = i;
        }                
    }
    return index;
}
int min_volume(struct solution_struct* solution){
    int cur = 0, min = 0;
    int index;
    for (int d = 0; d < solution->prob->dim; d++){
        min += solution->prob->items[0].size[d];  //initialize min
    }
    //printf("init min: %d\n", min);
    
    for (int i = 1; i < solution->prob->n; i++){
        for (int d = 0; d < solution->prob->dim; d++){
            cur += solution->prob->items[i].size[d];
        } 
        //printf("cur: %d\n", cur);
        if (min > cur && solution->x[i]== 0){
            min = cur;
            index = i;
            //printf("index %d, min: %d\n", index, min);
        }   
        cur = 0;             
    }
    return index;
}
//modify the solutions that violate the capacity constraints
void feasibility_repair(struct solution_struct* pop)
{
    //drop
    for (int p = 0; p < POP_SIZE; p++){
        evaluate_solution(&pop[p]);
        while (pop[p].feasibility < 0){
            for (int n = 0; n < pop[p].prob->n; n++){
                if(pop[p].x[n] == 1){
                    pop[p].x[n] = 0;
                    //printf("take out: %d\n", n);
                    break;
                }
            }
            evaluate_solution(&pop[p]);
        }   
        //printf("feas? %d\n", pop[p].feasibility);
    }
    //add
    
    /**
    for (int p = 0; p < POP_SIZE; p++){
        evaluate_solution(&pop[p]);
        while (pop[p].feasibility >= 0)
        {
            int minv = min_volume(&pop[p]);
            printf("min %d\n", minv);
            pop[p].x[minv] = 1;
            evaluate_solution(&pop[p]);
            //printf("feas? %d\n", pop[p].feasibility);
            if (pop[p].feasibility < 0)
            {
                pop[p].x[minv] = 0;
                //printf("feas: %d\n", pop[p].feasibility);
                break;
            }
        }
    }**/
    
}

//local search
void local_search_first_descent(struct solution_struct* pop)
{
    bool canPut = false, finished = false;
    for (int p = 0; p < POP_SIZE; p++){   //for each solution 
        for (int i = 0; i < pop[p].prob->n; i++){   //for each item 
            if (pop[p].x[i] == 1){   //if the item is in the knapsack
                for (int j = 0; j < pop[p].prob->n; j++){   
                    //find another item other than i, not in the sack
                    if (i != j && pop[p].x[j] == 0){
                        for (int d = 0; d < pop->prob->dim; d++){ 
                            //can fit in the sack in all dimensions and has more profit, i.e. a better neighbour
                            if(pop[p].cap_left[d] + pop[p].prob->items[i].size[d] > pop[p].prob->items[j].size[d] 
                            && pop[p].prob->items[i].p < pop[p].prob->items[j].p){
                                canPut = true;
                            }else{
                                canPut = false;
                            }   
                        }
                    }
                    if (canPut){
                        pop[p].x[i] = 0;
                        pop[p].x[j] = 1;
                        //printf("population: %d, out item: %d, in item: %d\n", p, i, j);
                        finished = true;
                        break;
                    }
                }
                if (finished)
                    break;
            }
        }
    }
}
//replacement
void replacement(struct solution_struct* curt_pop, struct solution_struct* new_pop)
{
    for (int i = 0; i < POP_SIZE; i++)
    {
        copy_solution(&curt_pop[i], &new_pop[i]);
    }
}

//update global best solution with best solution from pop if better
void update_best_solution(struct solution_struct* pop)
{
    int bestFit = 0, bestFitIndex = 0;
    for (int p = 0; p < POP_SIZE; p++)
    {
        evaluate_solution(&pop[p]);
        if (pop[p].objective > bestFit){
            bestFitIndex = p;
            bestFit = pop[p].objective;
            printf("best index = %d, best = %d\n", bestFitIndex, bestFit);
        }
    }
    copy_solution(best_sln, &pop[bestFitIndex]);
}
//memetic algorithm
int MA(struct problem_struct* prob)
{
    struct solution_struct curt_pop[POP_SIZE];
    struct solution_struct new_pop[POP_SIZE];
    init_population(prob, curt_pop);
    init_population(prob, new_pop);

    int gen=0;
    clock_t time_start, time_fin;
    time_start = clock();
    double time_spent=0;
    while(gen < MAX_NUM_OF_GEN && time_spent < MAX_TIME)
    {
        selection(curt_pop);
        cross_over(new_pop, curt_pop);
        mutation(new_pop);

        feasibility_repair(new_pop);
        local_search_first_descent(new_pop);
        replacement(curt_pop, new_pop);
        gen++;
        time_fin=clock();
        time_spent = (double)(time_fin - time_start)/CLOCKS_PER_SEC;
    }
    
    update_best_solution(curt_pop);
    
    free_population(curt_pop, POP_SIZE);
    free_population(new_pop, POP_SIZE);
    
    return 0;
}

int main(int argc, const char * argv[]) {
    
    printf("Starting the test run!\n");
    if (argc < 3){
        printf("too few arguments!\n");
        return -1;
    }
    if (argc > 9){
        printf("too many arguments!\n");
        return -2;
    }

    char out_file[20];
    strcpy(out_file, argv[4]);
    MAX_TIME = atoi(argv[6]);
    
    int num_of_problems;
    FILE* fp = fopen(argv[2], "r");
    FILE* fp1 = fopen(argv[4], "a");
    fscanf(fp, "%d", &num_of_problems);
    struct problem_struct** my_problems = load_problems(num_of_problems, fp);

    //allocate memory for best_sln
    best_sln = malloc(sizeof(struct solution_struct));
    best_sln->prob = malloc(sizeof(struct problem_struct));
    best_sln->cap_left = malloc(sizeof(int)*5);
    best_sln->x = malloc(sizeof(int)*POP_SIZE);
    
    for(int k=0; k<num_of_problems; k++)
    {
        for(int run=0; run<NUM_OF_RUNS; run++)
        {
            srand(RAND_SEED[run]);
            MA(my_problems[k]); //call MA
            output_solution(best_sln, out_file);
        }
        if (my_problems[k] != NULL)
        {
            free_problem(my_problems[k]); //free problem data memory
        }
        
    }
    if (my_problems != NULL)
    {
        free(my_problems); //free problems array
    }
    
    if(best_sln->x!=NULL && best_sln->cap_left!=NULL)
    { 
        free(best_sln->cap_left); 
        free(best_sln->x);
        if (best_sln!=NULL)
        {
            free(best_sln);
        }
    }
    
    fclose(fp);
    return 0;
}