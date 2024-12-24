#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <math.h>

#define EXIT 0
#define INVALID -1
#define CPU 1
#define PID 5
#define RACE 2
#define BANKER 3
#define PROCESS 5
#define	RESOURCE 3
#define TRUE 1
#define FALSE 0
#define MEMORY 4
#define FIRST 0
#define BEST 1
#define WORST 2
#define NEXT 3
#define PAGE 5
#define FIFO 0
#define LRU 1
#define FRAMES 4
#define DISK 6
#define FCFS 0
#define SSTF 1
#define REQUEST 8

int displayMenu(); //display options
void cpuScheduling(); //data input and arithmetic function calls
void fcfs(int * process, int * arrvTime, int * burstTime); // calculates values and calls the display functions
void sjf(int * process, int * arrvTime, int * burstTime); // similar function as fcfs with different scheduling calculation
void displaySchedule(int * process, int * at, int * bt, int * wt, int * tat); //display the numbers calculated
void raceCondition();
void * threadFuncOne();
void * threadFuncTwo();
void bankersAlgorithm();
void systemState(int * feasible, int * safe);

void memoryManagement();
void displayProcess(int * allocation, int processes, int * processSize);
void firstFit(int * blockSize, int blocks, int * processSize, int processes);
void worstFit(int * blockSize, int blocks, int * processSize, int processes);
void bestFit(int * blockSize, int blocks, int * processSize, int processes);
void nextFit(int * blockSize, int blocks, int * processSize, int processes); // TODO

void pageReplacement(); // page replacement algorithms
void fifo();
void lru();
int findLRU(int * time);
void displayPages();

void diskScheduling();
void diskFcfs(int * requests, int head);
void diskSstf(int * requests, int head);

//global variables
int resource = 5;

int main(){ 
	int choice = INVALID;
	while (choice != EXIT){ // while not exitting (0 is not inputted)
		choice = displayMenu();
		if (choice == CPU) cpuScheduling(); // if 1 is inputted, calculate and display the 2 scheduling algorithms
		else if (choice == EXIT) exit(choice); //exit if 0 is in play (ends loop)
		else if (choice == RACE) raceCondition(); // if 2 is inputted, display race condition simulation
		else if (choice == BANKER) bankersAlgorithm(); // if 3 is inputted, display bankersAlgorithm
		else if (choice == MEMORY) memoryManagement();
		else if (choice == PAGE) pageReplacement(); // if 5 is inputted, display the page replacement algorithms
		else if (choice == DISK) diskScheduling();
	}
	return 0;
}

int displayMenu(){
	int choice = INVALID;
	while (choice == INVALID){ // allows for the menu to be displayed initially
		printf("\n\n\n\t\t\t\t\t ********** Operating System Management Menu ********** \n");
		printf("\nSelect the OS program to run, enter the number of your selection.\n1. CPU Scheduling \n2. Race Condition\n3. Banker's Algorithm \n4. Memory Management\n5. Page Replacement \n6. Disk Scheduling\n0. Exit\n");
		// it was unclear to me if you wanted it to loop until 0 was inputted, but I set it up so it would be 
		scanf("%d", &choice); 
		if (choice > DISK || choice < EXIT) choice = INVALID; // if the value is outside the bounds of 0-3 (invlaid) FIXME when we add more inputs in future assignments
		}
	return choice;
}

void cpuScheduling(){ 
	int process[PID] = {1, 2, 3, 4, 5};
	int arrvTime[PID] = {0, 2, 4, 6, 7};
	int burstTime[PID] = {8, 5, 10, 2, 3};

	fcfs(process, arrvTime, burstTime); // call functions
	sjf(process, arrvTime, burstTime);	
}

void fcfs(int * process, int * at, int * bt){ // first come first serve algorithm
	int wt[PID] = {0};
	int tat[PID] = {0};
	wt[0] = 0;
	for (int i = 1; i <= PID; i++){ //calculate wait time by adding burst time and the wait time of the last element
		wt[i] = bt[i-1] + wt[i-1];
	}

	for (int i = 1; i <= PID; i++){ // calculate turn around time by adding burst time and waiting time of the last elements
		tat[i-1] = bt[i-1] + wt[i-1];
	}
	printf("\tFCFS\n");
	displaySchedule(process, at, bt, wt, tat); // display results
}

void sjf(int * process, int * at, int * bt){ // shortest job first algorithm
	int wt[PID] = {0};
	int tat[PID] = {0};
	int idx, temp;
	//sorting loop -- need to sort loop first since the algorithm is reliant on who has the fastest burst time
	for (int i = 0; i < PID; i++){ 
		idx = i;
		for (int j = i+1; j < PID; j++){
			if (bt[j] < bt[idx]) idx = j; // if element in burst time at spot idx is larger than at spot j, replace idx's value (max value) with current value (j)
		}
		temp = bt[i];
		bt[i] = bt[idx];
		bt[idx] = temp;

		temp = process[i];
		process[i] = process[idx];
		process[idx] = temp;

		temp = at[i];
		at[i] = at[idx];
		at[idx] = temp; // swapping variables to be held at location idx rather than i
	}
	wt[0] = 0;
	for (int i = 1; i <= PID; i++){ // calculating wait time
		wt[i] = bt[i-1] + wt[i-1];
	}
	for (int i = 1; i <= PID; i++){ // calculating burst time
		tat[i-1] = bt[i-1] + wt[i-1];
	}
	printf("\tSJF\n");
	displaySchedule(process, at, bt, wt, tat); // display results
}

void displaySchedule(int * process, int * at, int * bt, int * wt, int * tat){
	int totWt = 0;
	int totTat = 0;
	float avgWt = 0;
	float avgTat = 0;
	printf("PID AT BT WT TAT\n");
	for (int i = 0; i < PID; i++){ // print each row of data
		totWt += wt[i]; // sum of waiting times and turn around times
		totTat +=tat[i];
		printf("%2d  %2d %2d %2d %2d\n", process[i], at[i], bt[i], wt[i], tat[i]);
	}
	
	avgWt = (float)totWt / PID; // find averages and print
	avgTat =(float)totTat / PID;
	printf("Average waiting time = %.2f \nAverage turn around time = %.2f \n", avgWt, avgTat);
}	

void raceCondition(){
	pthread_t threadOne;
	pthread_t threadTwo;
	
	pthread_create(&threadOne, NULL, threadFuncOne, NULL);
	pthread_create(&threadTwo, NULL, threadFuncTwo, NULL);
	
	pthread_join (threadOne, NULL);
	pthread_join (threadTwo, NULL);

	printf("The value of the shared resource is: %d \n", resource);	
}

void * threadFuncOne(){
	int threadId = 1;
	int fOne = resource;

	printf("Thread 1 reads in value as: %d \n", fOne);
	fOne +=1;
	printf("Local update to shared resource updates the value to: %d \n\n", fOne);
	
	sleep(1);
	resource = fOne;

	printf("The value of shared resource by Thread 1 is: %d \n", resource);
}

void * threadFuncTwo(){
	int threadId = 2;
	int fTwo = resource;

	printf("Thread 2 reads in value as: %d \n", fTwo);
	fTwo -=1;
	printf("Local update to shared resource updates the value to: %d \n\n", fTwo);
	
	sleep(1);
	resource = fTwo;

	printf("The value of shared resource by Thread 2 is: %d \n", resource);
}

void bankersAlgorithm(){
	//variable and array initialization/setup
	int allocation[PROCESS][RESOURCE] = {{0, 0, 2},
										 {3, 0, 2},
										 {0, 1, 0},
										 {2, 1, 1},
									 	 {2, 0, 0}}; // 2d array for allocation
	
	int maxDemand[PROCESS][RESOURCE] = {{4, 3, 3},
										{9, 0, 2},
										{7, 5, 3},
										{2, 2, 2},
										{3, 2, 2}}; //  2d array for maximum resource demand

	int available[RESOURCE] = {2, 4, 6}; //  array for available resources

	int need[PROCESS][RESOURCE]; // array for needed allocations calculation
	int feasible[PROCESS]; // array for storing feasible resource allocation

	int safe[PROCESS];
	int safeIdx = 0; // index variable for traversing safe array
	
	for (int p = 0; p < PROCESS; p++){
		feasible[p] = FALSE; // initializing all values in feasible to false
	}

	for(int p = 0; p < PROCESS; p++){
		for(int r = 0; r < RESOURCE; r++){
				need[p][r] = maxDemand[p][r] - allocation[p][r]; // calculate needed resources using maxDemand-allocation arrays
		}
	}

	// BANKERS ALGORITHM
	for (int i = 0; i < PROCESS; i++){
		for (int p = 0; p < PROCESS; p++){
			if (feasible[p] == FALSE){
				int isUnsafe = FALSE;
				for (int r = 0; r < RESOURCE; r++){
					if (need[p][r] > available[r]){
						isUnsafe = TRUE;
						break;
						}
					}
				
				if (isUnsafe == FALSE){
					safe[safeIdx++] = p; // save the safe state for the process
					for (int r = 0; r < RESOURCE; r++){
						available[r] = available[r] + allocation[p][r];
					}
					feasible[p] = TRUE;
				}
			}
		}
	}
	systemState(feasible, safe);
}

void systemState(int * feasible, int * safe){
	int isSafe = TRUE;
	for (int i = 0; i < PROCESS; i++){
		if(feasible[i] == FALSE){
			isSafe = FALSE;
			printf("\n\nThe operating system state is unsafe");
			break;
		}
	}
	if (isSafe == TRUE){
		printf("\n\n");
		for (int i = 0; i < PROCESS; i++){
			printf("P%d", safe[i]);
			if (i < PROCESS-1) printf("->");
		}
	}
}

void memoryManagement(){
    int algorithm = FIRST;
	
    while (algorithm != NEXT+1) { // loop through each algorithm first to next (0-3)
		int blockSize[] = {70, 20, 45, 65, 40, 80};
		int processSize[] = {15, 35, 25, 45, 60, 20};
    	int blocks = sizeof(blockSize) / sizeof(blockSize[0]); // defines how many blocks there are in blockSize
    	int processes = sizeof(processSize) / sizeof(processSize[0]); // defines how many blocks there are in processSize

		if (algorithm == FIRST){
            firstFit(blockSize, blocks, processSize, processes);
		}
		else if(algorithm == BEST){
			bestFit(blockSize, blocks, processSize, processes);
        }
		else if (algorithm == WORST){
			worstFit(blockSize, blocks, processSize, processes);
		}
        else if (algorithm == NEXT){
			nextFit(blockSize, blocks, processSize, processes);
       	}
		algorithm++;
    } 
}
void nextFit(int blockSize[], int blocks, int processSize[], int processes) {
    int allocation[processes];
	int id = 0; // Variable to store the block allocation for a process, initialized to 0

    memset(allocation, INVALID, sizeof(allocation)); // Initialize allocation array with -1 (INVALID)

    for (int i = 0; i < processes; i++) {
    	//int id = 0;         
		while (id < blocks) {
            if (blockSize[id] >= processSize[i]) {
                allocation[i] = id; // Update the allocation array
                blockSize[id] -= processSize[i]; // Reduce available memory of the current block
                break; // Break out of the inner loop
            }
			id = (id + 1) % blocks; // Move to the next block	
        }
    }
	printf("\n\n\t\tNext Fit\n");
    displayProcess(allocation, processes, processSize); // Display process allocation
}

void firstFit(int blockSize[], int blocks, int processSize[], int processes) {
    int allocation[processes];

    memset(allocation, INVALID, sizeof(allocation)); // Initialize allocation array with -1 (INVALID)

    for (int i = 0; i < processes; i++) {
        for (int j = 0; j < blocks; j++) {
            if (blockSize[j] >= processSize[i]) {
                allocation[i] = j; // Update the allocation array
                blockSize[j] -= processSize[i]; // Reduce available memory of the current block
                break; // Break out of the inner loop
            }
        }
		
    }
    printf("\n\n\t\tFirst Fit\n");
	displayProcess(allocation, processes, processSize); // Display process allocation
}

void bestFit(int blockSize[], int blocks, int processSize[], int processes) {
    int allocation[processes];

    memset(allocation, INVALID, sizeof(allocation)); // Initialize allocation array with -1 (INVALID)

    for (int i = 0; i < processes; i++) {
        int bestIdx = INVALID; // Variable to store the current best fit value

        for (int j = 0; j < blocks; j++) {
            if (blockSize[j] >= processSize[i]) {
                if (bestIdx == INVALID || blockSize[j] < blockSize[bestIdx]) {
                    bestIdx = j; // Update the bestIdx if the current block is a better fit
                }
            }
        }
        if (bestIdx != INVALID) {
            allocation[i] = bestIdx; // Update the allocation array
            blockSize[bestIdx] -= processSize[i]; // Reduce available memory of the current block
        }
    }
	printf("\n\n\t\tBest Fit\n");
    displayProcess(allocation, processes, processSize); // Display process allocation
}

void worstFit(int blockSize[], int blocks, int processSize[], int processes) {
    int allocation[processes];

    memset(allocation, INVALID, sizeof(allocation)); // Initialize allocation array with -1 (INVALID)

    for (int i = 0; i < processes; i++) {
        int wstIdx = INVALID; // Variable to store the current worst fit value

        for (int j = 0; j < blocks; j++) {
            if (blockSize[j] >= processSize[i]) {
                if (wstIdx == INVALID || blockSize[j] > blockSize[wstIdx]) {
                    wstIdx = j; // Update the wstIdx if the current block is a worse fit
                }
            }
        }
        if (wstIdx != INVALID) {
            allocation[i] = wstIdx; // Update the allocation array
            blockSize[wstIdx] -= processSize[i]; // Reduce available memory of the current block
        }
    }
	printf("\n\n\t\tWorst Fit\n");
    displayProcess(allocation, processes, processSize); // Display process allocation
}

void displayProcess(int * allocation, int processes, int * processSize) {
    printf("\nProcess No.\tProcess Size\tBlock No.\n");
    for (int i = 0; i < processes; i++) {
        printf("%d\t\t%d\t\t", i+1, processSize[i]);
        if (allocation[i]+1 != INVALID+1) {
			printf("%d\n", allocation[i]+1);
        } else {
			printf("Not Allocated\n");
        }
    }
}


void pageReplacement(){
	for (int algorithm = 0; algorithm < 2; algorithm++) { // call each page replacement algorithm
		if (algorithm == FIFO)
			fifo();
		if (algorithm == LRU)
			lru();
	}
}

void fifo(){
	printf("\n******** First In First Out ********\n");
	printf("Page\tFrame 1\tFrame 2\tFrame 3\t Frame 4\t\n");
	int pageRequests[15] = {2, 3, 8, 4, 5, 6, 5, 7, 1, 8, 3, 1, 4, 2, 6};
	int pageFaults = 0;
	int allocation[FRAMES];
	int present;
	int pages = sizeof(pageRequests) / sizeof(pageRequests[0]);	
	memset(allocation, INVALID, sizeof(allocation));

	for (int i = 0; i < pages; i++){ // for the number of pages
		present = 0;
		for (int j = 0; j < FRAMES; j++){ // ...and for each frame
			if (pageRequests[i] == allocation[j]){ // if the the request matches the allocaiton
				present++;
				pageFaults--;
			}
		}
		pageFaults++; // add a page fault if there is one
		if (pageFaults <= FRAMES && present == 0){ // if the there are more page faults then frames AND the present bit is still set to false
			allocation[i] = pageRequests[i]; //FIXED: changed allocation index to i
		}
		else if(present == 0){ // if just the present bit = false
			allocation[(pageFaults-1) % FRAMES] = pageRequests[i];
		}
		displayPages(pageRequests[i], allocation); // call display function
		printf("\n");
	}
	printf("\n");
	printf("Page Faults: %d \n", pageFaults); // print final page fault count
}

void lru(){
	printf("\n******** Least Recently Used ********\n");
	printf("Page\tFrame 1\tFrame 2\tFrame 3\t Frame 4\t\n");

	int pageRequests[15] = {2, 3, 8, 4, 5, 6, 5, 7, 1, 8, 3, 1, 4, 2, 6};
	int pageFaults = 0;
	int allocation[FRAMES];
	int pages = 15;
	int counter = 0;
	int time[10];
	int flag1;
	int flag2;
	int position = 0;
	memset(allocation, INVALID, sizeof(allocation));

	for(int i = 0; i < pages; i++){ // for number of pages
		flag1 = 0;
		flag2 = 0;
		for (int j = 0; j < FRAMES; j++){ // ...and number of frames in each page
			if(pageRequests[i] == allocation[j]){ 
				counter++;
				time[j] = counter;
				flag1 = 1;
				flag2 = 1;
				break;
			}
		}
		if (flag1 == 0){ // if the flag bit = 0
			for (int k = 0; k < FRAMES; k++){
				if (allocation[k] == INVALID){
					counter++;
					pageFaults++;
					allocation[k] = pageRequests[i];
					time[k] = counter;
					flag2 = 1;
					break;
				}
			}
		}
		if (flag2 == 0){ // if the flag bit is 0:
			position = findLRU(time);
			counter++;
			pageFaults++;
			allocation[position] = pageRequests[i]; // set the allocation = to pagerequest
			time[position] = counter; // save the counter variable to time at position 
		}
		displayPages(pageRequests[i], allocation);
		printf("\n");
	}
	printf("\n"); //FIXED: added new line for formatting?
	printf("Page Faults: %d \n", pageFaults);
}

int findLRU(int * time){ // simple function to find the last recently used frame
	int position = 0;
	int minimum = time[0];

	for (int i = 0; i < FRAMES; i++){
		if(time[i] < minimum){
			minimum = time[i];
			position = i;
		}
	}
	return position; // return where it was
}

void displayPages(int page, int * allocation){ // simple display function to show page data of allocations
	printf("%d ", page);
	for (int i = 0; i < FRAMES; i++){
		if (allocation[i] == INVALID) printf("-");
		else printf("\t%d ", allocation[i]);
	}
}

void diskScheduling(){
	int requests[REQUEST] = {146, 89, 24, 70, 102, 13, 51, 134}; // stores disk requests
	int head = 50; // where the disk location counter starts
	int algorithm = 0;
	printf(" ******** Disk Scheduling  ********\n\n\n");
	while (algorithm < SSTF+1){
		if(algorithm == FCFS)
			diskFcfs(requests, head);
		else if (algorithm == SSTF)
			diskSstf(requests, head);
		algorithm++;
	}
}

void diskFcfs(int * requests, int head){
	int seek = 0; // number of times a seek operation was performed
	int track = 0; // current location of the track
	int distance = 0; // total seek distance
	int start = head; // stores where the head starts, thus initialized to head


	printf("  ******** First Come First Serve  ******** \n\n");
	printf("Seek Sequence: %d", head);
	for (int i=0; i < REQUEST; i++){ // for each request...
		track = requests[i];
		distance = abs(head - track); // FIXME could have issues
		seek += distance;
		head = track;
	}
	for (int i=0; i < REQUEST; i++){ // for each request...
		printf(" -> %d", requests[i]);
	}
	printf("\nTotal Seek Operations: %d\n\n", seek);
}

void diskSstf(int * requests, int head){
	int sequence[REQUEST] = {0}; // store the sequence of service requests
	int distance[REQUEST] = {0}; // store distance between service requests
	int seek = 0; // store the number of seek operations
	int start = head; // the starting location of the "needle"
	int minVal = 0; // the minimum distance value
	int minValIdx = 0; // the index of the minimum value
	int seqIdx = 0; // store the index of the sequence


	for (int i=0; i < REQUEST; i++){ // for each request...
		for (int j=0; j < REQUEST; j++){ // for each request... calculate distance between each request
			distance[j] = abs(head - requests[j]); // store the value of distance for each request
		}
		minVal = distance[0]; // set the mininum value to the arbitrary default
		minValIdx = 0; // set the minimum value to 0;

		for (int k = 1; k < REQUEST; k++){ // find minimum distance
			if (minVal > distance[k]){
				minVal = distance[k];
				minValIdx = k;
			}
		}
		sequence[seqIdx++] = requests[minValIdx]; // NOTE: not what was asked for (?) but I got it to work this way
		head = requests[minValIdx]; // set the value head equal to the next indx in requests
		requests[minValIdx] = 999; // don't process this minVal again: set to very large value
	}
	printf("  ******** Shortest Seek Time First ******** \n\n");
	printf("Seek Sequence: %d ->", start);
	seek += abs(start-sequence[0]);
	printf(" %d", sequence[0]);

	for (int l = 1; l < REQUEST; l++){
		seek += abs(sequence[l]-sequence[l-1]);
		printf(" -> %d", sequence[l]);
	}
	printf("\nTotal Seek Operations: %d", seek);
}

