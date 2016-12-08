//Coarse and Fine Grained - Conditional Compilation
//Enter no of threads at runtime
//g++ -std=c++11 -DUSE_COARSE locks.cpp -lpthread
//g++ -std=c++11 -DUSE_FINE locks.cpp -lpthread

#include <stdio.h>
#include <iostream>
#include <pthread.h>
#include <stdlib.h>
#include <list>
#include <vector>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <random>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>

using namespace std;

#define DIGITS 4

fstream myfile("barrier_live_journal.txt", std::ios_base::in);
ofstream outfile("output.txt",std::ios_base::app);

int numver, max_deg=0, p;

#ifdef USE_COARSE
pthread_mutex_t lock;
#endif

typedef struct 
{
    int     secs;
    int     usecs;
}TIME_DIFF;

struct node
{
	int id;
	node *next;
};

struct Vertex
{
	bool boundary_vertex;
	int Vid;
	int Partition_id;
	int color;
	node *next;
	#ifdef USE_FINE
	pthread_mutex_t lock;
	#endif
};

Vertex *ver;

TIME_DIFF * my_difftime (struct timeval * start, struct timeval * end)
{
	TIME_DIFF * diff = (TIME_DIFF *) malloc ( sizeof (TIME_DIFF) );
 
	if (start->tv_sec == end->tv_sec) 
	{
        	diff->secs = 0;
        	diff->usecs = end->tv_usec - start->tv_usec;
    	}
   	else 
	{
        	diff->usecs = 1000000 - start->tv_usec;
        	diff->secs = end->tv_sec - (start->tv_sec + 1);
        	diff->usecs += end->tv_usec;
        	if (diff->usecs >= 1000000) 
		{
        	    diff->usecs -= 1000000;
	            diff->secs += 1;
	        }
	}
        return diff;
}

void input() //Input should be edges - undirected or directed
{
	int d,numedge,temp1,temp2,edge=0;
	myfile>>numedge;
	node *ptr;
	for(int i=1;i<=numver;i++)
	{
		ver[i].Partition_id = -1;
		ver[i].next = NULL;
		ver[i].Vid = i;
		ver[i].boundary_vertex = false;
		ver[i].color = -1;
	}

	for(int i=0;i<numedge;i++)
	{
//		cout<<i<<endl;
		myfile >> temp1;
		myfile >> temp2;
		if(ver[temp1].next == NULL)
		{
			ver[temp1].next=(node*)(malloc(sizeof(node)));
			ptr = ver[temp1].next;
			ptr->id = temp2;
			edge++;
			ptr->next = NULL;
		}
		else
		{
			ptr = ver[temp1].next;
			while(ptr->next != NULL)
			{
				if(ptr->id == temp2)
					break;			
				ptr = ptr->next;
			}
			if(ptr->next == NULL && ptr->id != temp2)
			{
				ptr->next = (node*)(malloc(sizeof(node)));
				ptr = ptr->next;
				ptr->id = temp2;
				edge++;				
				ptr->next = NULL;
			}	
		}
  	        if(ver[temp2].next == NULL)
		{
			ver[temp2].next=(node*)(malloc(sizeof(node)));
			ptr = ver[temp2].next;
			ptr->id = temp1;
			ptr->next = NULL;
		}
		else
		{
			ptr = ver[temp2].next;
			while(ptr->next != NULL)
			{
				if(ptr->id == temp1)
					break;			
				ptr = ptr->next;
			}
			if(ptr->next == NULL && ptr->id != temp1)
			{
				ptr->next = (node*)(malloc(sizeof(node)));
				ptr = ptr->next;
				ptr->id = temp1;
				ptr->next = NULL;
			}
		}
	}
	cout<<"No of undirected edges: "<< edge << endl;
	outfile<<"No of undirected edges: "<< edge<< endl;
	for(int i=1; i<=numver; i++) 
	{
		d=0;
		ptr = ver[i].next;
		while(ptr != NULL)
		{
			d++;
			ptr=ptr->next;
		}		
		if(d > max_deg)
			max_deg = d;
	}
	outfile<<"Max Degree: "<< max_deg << endl;
}

void* color(void* t) 
{
	long tid=(long)t;

	list<Vertex*> ulist; //all internal vertices of this thread
	list<Vertex*> blist; //boundary vertices of this thread

	int k;
	int avclr[max_deg+1];
	list<int> temp,t1;

	for(int j=0; j < max_deg+1; j++) 
		avclr[j]=j;

	long partition_id;

	for(int i=1;i<=numver;i++)
	{
		partition_id = (long)ver[i].Partition_id;
//		cout<<partition_id<<" "<<tid<<" "<<ver[i].Vid<<" "<<ver[i].boundary_vertex<<endl;
		if(partition_id == tid)
		{
			if(ver[i].boundary_vertex == 0)
				ulist.push_back(&ver[i]);
			else
				blist.push_back(&ver[i]);
		}
	}

	int id,temp1;
	node* tp;

	for(list<Vertex*>::iterator iter=ulist.begin(); iter!=ulist.end(); iter++) 
	{
		temp.clear();
		tp=(*iter)->next;

		while(tp!=NULL) //Check the colors of adjacent vertices and assign a least possible color different from the adjacent vertices.
		{
			id=tp->id;
			if(ver[id].boundary_vertex == false) 
			{
				if(ver[id].color!=-1)
					temp.push_back(ver[id].color);
			}
			tp=tp->next;
		}
		k=0;
		list<int>::iterator t = temp.begin();
		while(t != temp.end())
		{
			if(*t == k)
			{
				k++;
				t = temp.begin();
				continue;
			}
			t++;
		}
		(*iter)->color = k;
	}
	cout<<tid<<"'s Internal Vertices Colored"<<endl;

	for(list<Vertex*>::iterator iter=blist.begin(); iter!=blist.end(); iter++) 
	{
		temp.clear();

		#ifdef USE_COARSE
		pthread_mutex_lock(&lock);
		#endif
		
		#ifdef USE_FINE
		t1.clear();
		t1.push_back((*iter)->Vid);
		tp = (*iter)->next;
		while(tp!=NULL)
		{
			temp1 = tp->id;
			if(ver[temp1].boundary_vertex == true)
				t1.push_back(temp1);  //to lock all adjacent boundary vertices
			tp = tp->next;		
		}		
		t1.sort();
		
		for(list<int>::iterator t=t1.begin(); t!=t1.end(); t++)
			pthread_mutex_lock(&ver[*t].lock);
		#endif

		tp=(*iter)->next;

		//Check the colors of adjacent vertices and assign a least possible color different from the adjacent vertices.
		while(tp!=NULL) 
		{
			id=tp->id;
			if(ver[id].color!=-1)
				temp.push_back(ver[id].color);
			tp=tp->next;
		}
		k=0;
		list<int>::iterator t = temp.begin();
		while(t != temp.end())
		{
			if(*t == k)
			{
				k++;
				t = temp.begin();
				continue;
			}
			t++;
		}
		(*iter)->color = k; //color the boundary vertex

		#ifdef USE_COARSE
		pthread_mutex_unlock(&lock);
		#endif

		#ifdef USE_FINE
		for(list<int>::iterator t=t1.begin(); t!=t1.end(); t++)
			pthread_mutex_unlock(&ver[*t].lock);
		#endif
	}
	cout<<tid<<"'s boundary vertices completed"<<endl;
}

int main(int argc,char *argv[]) 
{
	myfile >> numver;
	cout<<"The no of vertices: "<<numver<<endl;

	#ifdef USE_COARSE	
	outfile<<"COARSE GRAINED LOCKING"<<endl;
	#elif USE_FINE	
	outfile<<"FINE GRAINED LOCKING"<<endl;
	#else
	cout<<"Error: Use a directive (Coarse/Fine) to compile"<<endl;	
	return 0;
	#endif
		
	outfile<<"The no of vertices: "<<numver<<endl;
	int NTHREADS;

	if(argc > 1)
		NTHREADS = atoi(argv[1]);
	else
	{
		cout<<"Specify No of Threads: ";
		cin >> NTHREADS;
	}

	if(numver < NTHREADS)
	{
		cout<<"Error: No of threads is more than no of vertices"<<endl;
		return 0;
	}	

	ver = (Vertex*) malloc((numver+1) * sizeof(Vertex));

	input();
		
	#ifdef USE_COARSE
	pthread_mutex_init(&lock,NULL);

	#elif USE_FINE
	for(int i=1;i<=numver;i++)
		pthread_mutex_init(&ver[i].lock,NULL);
	#endif

	double duration;
 	time_t start1, start2;
	time_t end1, end2;
   
	int i,j,no_colors=0;
	int temp,dig;
	node* pt;
	double b;
	
	cout<<"Graph Generated"<<endl;
	//Printing the graph
/*
	node *ptr;
	for(int i=1; i<=numver; i++) 
	{
		cout<<ver[i].Vid<<": ";
		ptr=ver[i].next;
		while(ptr!=NULL) 
		{
			cout<<ptr->id<<" ";
			ptr=ptr->next;
		}
		cout<<"Partition: "<<ver[i].Partition_id<<" ";    //Partition ids are from 0,....,NTHREADS-1
		cout<<"Boundary vertex: "<<ver[i].boundary_vertex<<" ";
		cout<<"Color: "<<ver[i].color;
		cout<<endl;
	}
	cout<<"Max Degree: "<<max_deg<<endl;
*/

	pthread_t *thr = new pthread_t[NTHREADS];
    	// Make threads Joinable for sure.
    	pthread_attr_t attr;
   	pthread_attr_init (&attr);
   	pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
	struct timeval tv1, tv2;
	TIME_DIFF * difference;

//	start1 = clock();	
//	start2 = time (NULL);
//    	gettimeofday(&tv1, NULL);

	double a;
	vector<double> vec;

//	for(int th=1;th<=10;th++)  //average over 10 iterations
	{
		for(i=1;i<=numver;i++)
		{
			ver[i].color = -1;
			ver[i].Partition_id = -1;
			ver[i].boundary_vertex = false;
		}

		if(NTHREADS == 1)
		{
			for(i=1;i<=numver;i++)
			{
				ver[i].Partition_id = 0;
				ver[i].boundary_vertex = false;
			}
			b = 0.0;
		}

		if(NTHREADS > 1)
		{
			float x = (float)numver / (float)NTHREADS;
			p=floor(x); //no of vertices in 1 partition
	
			//Assigning Partition Id's randomly

			srand(unsigned(std::time(0)));
			vector<int> myvector;
			for (int i=1; i<=numver; ++i) 
				myvector.push_back(i);

//			vector<int>::iterator itemp=myvector.begin();
//			while(itemp != myvector.end())
//			{
//					cout<<*itemp<<" ";
//				itemp++;
//			}
//			cout<<"After shuffling";
	
			random_shuffle ( myvector.begin(), myvector.end() );
			vector<int>::iterator itemp=myvector.begin();
//			itemp = myvector.begin();
//			while(itemp != myvector.end())
//			{
//				cout<<*itemp<<" ";
//				itemp++;
//			}
//			cout<<"No of vertices in 1 thread: "<<p;
			int u=0;
	
 			struct timeval tv3, tv4;
			TIME_DIFF * difference1;

			gettimeofday(&tv3,NULL);
	
			for(i=1; i<numver; i=i+p)
			{
				j=i;
				while(j < (i+p) && j<=numver)
				{
					ver[*itemp].Partition_id = u;	
//					cout<<*itemp<<" "<<ver[*itemp].Partition_id<<endl;		
					j++;
					++itemp;
				}
				u++;
			}	
	
			for(i=1;i<=numver;i++)
				if(ver[i].Partition_id == -1)
					ver[i].Partition_id = rand()%NTHREADS;
	
			for(i=1; i<=numver; i++) 
			{
				pt=ver[i].next;
				while(pt!=NULL) 
				{
					temp=pt->id;
					if(ver[temp].Partition_id!=ver[i].Partition_id) 
					{
 		    				ver[i].boundary_vertex=true;
						break;
					}
					pt=pt->next;		
				}
			}
	
			gettimeofday(&tv4,NULL);
			difference1 = my_difftime (&tv3, &tv4);
			dig = 1;
			temp = difference1->usecs;
			b = 0.0;
	
			while(temp>=10)
			{	
				dig++;
				temp = temp/10;
			}
			temp =1;

			for(i=1;i<=dig;i++)
				temp = temp * 10;
			b = (double) difference1->secs + ((double)difference1->usecs / (double)temp);
		}


		gettimeofday(&tv1,NULL);
	
		for(int i=0; i<NTHREADS; i++) 
			pthread_create(&thr[i], &attr, color, (void*) i);

		for(int i=0; i<NTHREADS; i++) 
			pthread_join(thr[i],NULL);

		gettimeofday(&tv2,NULL);

		difference = my_difftime (&tv1, &tv2);
		dig = 1;
   		temp = difference->usecs;
	
		while(temp>=10)
		{	
			dig++;
			temp = temp/10;
		}
		temp =1;
		for(i=1;i<=dig;i++)
			temp = temp * 10;
		a = 0.0;
		a = (double) difference->secs + ((double)difference->usecs / (double)temp);
		a = a + b;		//adding time for partitioning

		vec.push_back(a);
	}

	double average = accumulate(vec.begin(),vec.end(),0.0) / vec.size();

 // 	gettimeofday(&tv2, NULL);
 //	difference = my_difftime (&tv1, &tv2);

	outfile<<"No of Threads : "<<NTHREADS<<endl;

//    	cout << "Duration 1 (clock function): " << fixed << setprecision (DIGITS) << static_cast<double> ((end1 - start1) / CLOCKS_PER_SEC) << " secs." << endl;
  //  	cout << "Duration 2 (time function): " << (end2 - start2) << " secs." << endl;
    //	cout << "Duration 3 (gettimeofday() function): " << difference->secs << "." << difference->usecs<<" secs."<<endl;

//	outfile<<"Time Taken : "<<(end2-start2)<<" secs."<<endl; 
    	cout << "Duration (gettimeofday() function) averaged over 10 runs: " << average <<" secs."<<endl;

    	outfile << "Duration (gettimeofday() function) averaged over 10 runs: " << average <<" secs."<<endl;

	for(i=1;i<=numver;i++)
	{
		if(ver[i].color > no_colors)
			no_colors = ver[i].color;
	}
	no_colors++;
	
//	for(i=1;i<=numver;i++)
//		cout<<ver[i].Vid<<" "<<ver[i].color<<endl;
	cout<<"No of colors used: "<<no_colors<<endl;
	outfile<<"No of colors used: "<<no_colors<<endl;

	free(ver);
	free(difference);
	free(thr);
	
	return 0;
}
