#include <cstdio> 
#include <cstdlib> 
#include <cilk/cilk.h> 
#include "timer.h" 
#include "io.h" 
#include <cilk/reducer_opadd.h> 

void initialize (int **World, int **gen,int N); // sets the arrays equal to SOMETHING, as of now it does a pattern 
//void printResult (int **gen,int N); // prints the array 
void calculate ( int **gen, int i, int j, int live); //calculates the NEXT generation 
void getNeighbor(int **World, int **gen,int N); //sees how many "neighbors" each "cell" has 
void reset_World (int **World, int **gen,int N); // resets World so it becomes equal to the NEXT
void initialize (int **World, int **gen,int N) 
{ 
    srand (time(NULL)); 
    for (int i = 0; i < N; i++) 
    { 
        for (int j= 0; j < N; j++) 
        { 
            World[i][j] = rand()%2; 
            gen[i][j] = rand()%2; 
        } 
    } 
} 
 
void initialize1(int **World, int **gen,int N) 
{ 
        for(int i = 0;i<N;i++) 
        { 
                for(int j=0;j<N;j++) 
                { 
                        gen[i][j]=World[i][j]; 
                } 
        } 
}       


void printResult(int **gen,int N) // printing the generation 
{ 
    for (int i = 0; i < N; i++) 
    { 
        for (int j = 0; j < N; j++) 
        { 
            cout.width(2); 
            cout << gen[i][j]; 
        } 
        cout << endl; 
    } 
} 


void getNeighbor(int **World, int **gen,int N) 
{ 
        int sum = 0; 
        cilk::reducer_opadd<int> live; 
        cilk_for (int i = 0; i<N;i++) { 
        for (int j = 0; j < N; j++) 
        { 
            if( ((i-1<N)&&(i-1>=0)&&(j-1<N)&&(j-1>=0))&&World[i-1][j-1] == 1) 
            { 
                live++; 
            } 
            if (((i<N)&&(i>=0)&&(j-1<N)&&(j-1>=0))&&World[i][j-1] == 1) 
            { 
                live++; 
            } 
            if (((i+1<N)&&(i+1>=0)&&(j-1<N)&&(j-1>=0))&&World[i+1][j-1] == 1) 
            { 
                live++; 
            } 
            if (((i-1<N)&&(i-1>=0)&&(j<N)&&(j>=0))&&World[i-1][j] == 1) 
            { 
                live++; 
            } 
            if(((i-1<N)&&(i-1>=0)&&(j+1<N)&&(j+1>=0))&&World[i-1][j+1] == 1) 
            { 
                live++; 
            } 
            if (((i<N)&&(i>=0)&&(j+1<N)&&(j+1>=0))&&World[i][j+1] == 1) 
            { 

		live++; 
            } 
            if (((i+1<N)&&(i+1>=0)&&(j+1<N)&&(j+1>=0))&&World[i+1][j+1] == 1) 
            { 
                live++; 
            } 
            if (((i+1<N)&&(i+1>=0)&&(j<N)&&(j>=0))&&World[i+1][j] == 1) 
            { 
                live++; 
            } 
            sum = live.get_value(); 
            live.set_value(0); 
            calculate(gen,i,j,sum); 
            //live = 0; 
        } 
         
    } 

} 


void calculate(int **gen, int i, int j, int live) 
{ 
         
    if((gen[i][j]== 1)&&(live <= 1)) //determine who "lives" "dies" and is "born" 
    { 
        gen[i][j] = 0; // dead 
    } 
         
    if ((gen[i][j]== 1) &&(live > 3 )) 
    { 
        gen[i][j] = 0; //die 
    } 
    if( (gen[i][j]== 1) && (live == 2 || 3)) 
    { 
        gen[i][j] = 1; //lives 
    } 
    if( (gen[i][j]== 0) && (live == 3)) 
    { 
        gen[i][j] = 1; 
    } 
} 

 

void reset_World(int **World, int **gen,int N) 
{ 
    cilk_for (int i = 0; i <N ; i++) 
    { 
        for (int j=0; j < N; j++) 
        { 
            World[i][j] = gen[i][j]; //sets World equal to the new generation 
        } 
    } 
} 




// Main method       
int main(int argc, char* argv[]) { 
        int N,M; 
        int **World; 
        double elapsedTime; 

        // checking parameters 
        if (argc != 3 && argc != 4) { 
                cout << "Parameters: <N> <M> [<file>]" << endl; 
                return 1; 
        } 
        N = atoi(argv[1]); 
        M = atoi(argv[2]); 
// allocating matrices 
        World = new int*[N]; 
        for (int i=0; i<N; i++){ 
                World[i] = new int[N]; 
        } 

        // reading files (optional) 
        if(argc == 4){ 
                //cout<<"You don't need to pass in the file. Try again!\n"; 
                //return 1; 
                readMatrixFile(World,N,argv[3]); 
        } 

        // starting timer 
        timerStart();


  // YOUR CODE GOES HERE 
        //if(argc == 3){ 
                int **gen; 
                gen = new int*[N]; 
                for (int i=0; i<N; i++){ 
                        gen[i] = new int[N]; 
                } 
                if(argc==3){ 
                        initialize (World,gen,N); 
                } 
                if(argc==4){ 
                        initialize1(World,gen,N); 
                } 
                int count =M; 
             
                while(M>0){ 
                        getNeighbor(World, gen,N); 
                        //printMatrix(gen,N); 
                        reset_World(World, gen,N); 
                        M = M -1; 
                        //cout << ""<<endl; 
                } 
        //} //testing the results is correct 
        if(argc == 4){ 
                printMatrix(World,N); 
        } 
        if(argc==3){ 
                printMatrix(World,N); 
        }        
        // stopping timer 
        elapsedTime = timerStop(); 

        cout << "Duration: " << elapsedTime << " seconds" << std::endl; 

        // releasing memory 
        for (int i=0; i<N; i++) { 
                delete [] World[i]; 
        } 
        delete [] World; 

        return 0;        
} 


