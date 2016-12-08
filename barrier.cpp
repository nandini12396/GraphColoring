//Barrier - Varying no of threads (Graph is fixed)
//Input graph must be from 1, ..., n

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
ofstream outfile("paper.txt",std::ios_base::app);

int numver, max_deg=0, p;

pthread_mutex_t lock;
pthread_barrier_t bar1,bar2;
int cnt = 0, NTHREADS;

typedef struct 
{
    int     secs;
    int     usecs;
}TIME_DIFF;

struct node 
{
	int id;
	int color;
	node *next;
};

struct Vertex 
{
	bool boundary_vertex;
	int Vid;
	int Partition_id;
	int color;
	node *next;
};

Vertex* ver;

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
			ptr->color = -1;
			ptr->next = NULL;
			edge++;
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
				ptr->color = -1;
				edge++;
				ptr->next = NULL;
			}	
		}
	        if(ver[temp2].next == NULL)
		{
			ver[temp2].next=(node*)(malloc(sizeof(node)));
			ptr = ver[temp2].next;
			ptr->id = temp1;
			ptr->color = -1;
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
				ptr->color = -1;
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

	list<Vertex*> ulist;     //all vertices of this thread
	list<Vertex*> recolor;

	int k;
	int avclr[max_deg + 1];

	for(int j=0; j < max_deg+1; j++)  //initialise available array
		avclr[j]=j;

	long partition_id;

	for(int i=1;i<=numver;i++)
	{
		partition_id = (long)ver[i].Partition_id;
		if(partition_id == tid)
			ulist.push_back(&ver[i]);
	}

	int temp1,ret,iteration = 0;
	list<int> t1;
	int i=0,size;
//	size = ulist.size();
	node *tp,*tp1;

	if(ulist.size() == 0)
	{
		cout<<"ERROR: "<< tid << " has ulist initially empty!" <<endl;
		exit(0);
	}

	while(true)          //make sure initially ulist has at least 1 vertex i.e. no partition is empty
	{
		pthread_mutex_lock(&lock);
		size = cnt;		
		pthread_mutex_unlock(&lock);

		if(size == NTHREADS)
			break;

//		if	ulist.size() != 0 || cnt != NTHREADS

		if(ulist.size() == 0)
			iteration = 1;

		for(list<Vertex*>::iterator iter=ulist.begin(); iter!=ulist.end(); iter++) 
		{
			tp = (*iter)->next;	
			t1.clear();
			while(tp != NULL)
			{
				if(tp->color != -1)
				t1.push_back(tp->color);
				tp = tp->next;		
			}			
			k=0;
			list<int>::iterator t = t1.begin();
			while(t != t1.end())
			{
				if(*t == k)
				{
					k++;
					t = t1.begin();
					continue;
				}
				t++;
			}
			(*iter)->color = k;
			tp = (*iter)->next;			//neighbours of internal vertices will be only boundary vertices of same partition
			while(tp!=NULL)
			{
				tp1 = ver[tp->id].next;
				
				while(tp1!=NULL)
				{
					if((tp1->id) == ((*iter)->Vid))
					{
						tp1->color = (*iter)->color;
						break;
					}
					tp1 = tp1->next;
				}
				tp = tp->next;				
			}
		}

		//Barrier to synchronise threads
		pthread_barrier_wait(&bar1);
//		cout<<tid <<"Reached here";

		recolor.clear();
		
		for(list<Vertex*>::iterator iter=ulist.begin(); iter!=ulist.end(); iter++) 
		{
			if((*iter)->boundary_vertex == true)	// boundary vertices
			{
				tp = (*iter)->next;             //neighbours of boundary vertex are across partitions
				while(tp!=NULL)
				{
				//	if(ver[tp->id].Partition_id != (*iter)->Partition_id)
					tp->color = ver[tp->id].color;
					tp = tp->next;
				}
/*

					tp1 = ver[tp->id].next;
					
					while(tp1!=NULL)
					{
						if((tp1->id) == ((*iter)->Vid))
						{
							tp1->color = (*iter)->color;
							break;
						}
						tp1 = tp1->next;
					}
					tp = tp->next;				
				}
*/				
				tp = (*iter)->next;
				while(tp!=NULL)
				{
					if(((*iter)->color == ver[tp->id].color) && ((*iter)->Partition_id < ver[tp->id].Partition_id))
						recolor.push_back(&ver[(*iter)->Vid]);
					tp = tp->next;
				}						
			}
		}
		cout << recolor.size() << endl;
		//copy recolor to ulist
		ulist.clear();
		ulist.assign( recolor.begin(),recolor.end() );

		if(ulist.size() == 0 && iteration == 0)	//check if colored all vertices of the partition
		{	
			cout << "round complete for 1 thread" << endl;
				pthread_mutex_lock(&lock);
				cnt++;
				pthread_mutex_unlock(&lock);
		}			
		//Barrier to synchronise threads
		pthread_barrier_wait(&bar2);
	//	if(recolor.size() == 0)
	//	cout << "Round complete."<<endl;
	}
}

int main(int argc,char *argv[]) 
{
	myfile>>numver;
	cout<<"The no of vertices: "<<numver<<endl;
	outfile<<"orkut DATA"<<endl;
	outfile<<"USING BARRIERS"<<endl;
	outfile<<"The no of vertices: "<<numver<<endl;

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

	ver=(Vertex*) (malloc((numver+1) * sizeof(Vertex)));
	input();

	int i,j,no_colors=0,id=0;
	pthread_mutex_init(&lock,NULL);

	pthread_barrier_init(&bar1,NULL,NTHREADS);
	pthread_barrier_init(&bar2,NULL,NTHREADS);

// 	time_t start1, start2;
//    	time_t end1, end2;
	
	int temp,dig;
	node* pt;

	pthread_t *thr = new pthread_t[NTHREADS];
    	// Make threads Joinable for sure.
    	pthread_attr_t attr;
   	pthread_attr_init (&attr);
   	pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
	struct timeval tv1, tv2;
	TIME_DIFF *difference;

	double a,b;
	vector<double> vec;

//	for(int th=1;th<=10;th++)  //average over 10 iterations
	{
		cnt = 0;

		for(i=1;i<=numver;i++)
		{
			ver[i].color = -1;
			ver[i].Partition_id = -1;
			ver[i].boundary_vertex = false;
			pt = ver[i].next;
			while(pt != NULL)
			{
				pt->color = -1;
				pt = pt->next;
			}
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
			p=(int)floor(numver/NTHREADS); //no of vertices in 1 partition

			//Assigning Partition Id's

			struct timeval tv3, tv4;
			TIME_DIFF * difference1;
		
			int u=0,z=0, i =1;

			gettimeofday(&tv3,NULL);
/*	
			for(int i=0; i<NTHREADS; i++) 
			{
				for(int j=0; ((j<p) && ((i*p+j) < numver)); j++) 
				{
					ver[i*p+j].Partition_id = i;
				}
			}
*/
//			ver[1].Partition_id = 0;
			
//			for(i=2; i<=numver; i = i+p)
			while(i <=  numver)
			{
//				z = u*i;
				for(j=0; j<p && i<=numver; j++)
//				while(j < (i+p) && j<=numver)
				{
					if(u == NTHREADS - 1)
					{
						while(i<=numver)
						{
							ver[i].Partition_id = u;		
							i++;
						}
					}
					ver[i].Partition_id = u;
					i++;
				}
				u++;
			}
	
			for(i=1;i<=numver;i++)
				if(ver[i].Partition_id == -1)
					ver[i].Partition_id = NTHREADS-1;

			for(i=1;i<=numver;i++)
				if(ver[i].Partition_id != ver[i+1].Partition_id)
				{
					cout << i << " " << ver[i].Partition_id << endl;
					cout << i+1 << " "<<ver[i+1].Partition_id << endl; 
				}

			//Identifying boundary vertices
	
//			int temp;
			cout << "partitioned";		
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
			cout << "boundary done";
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
			temp = 1;
	
			for(i=1;i<=dig;i++)
				temp = temp * 10;
			b = (double) difference1->secs + ((double)difference1->usecs / (double)temp);
		}

		gettimeofday(&tv1,NULL);
	
		for(int i=0; i<NTHREADS; i++) 
			pthread_create(&thr[i], &attr, color, (void*) i);

		for(int i=0; i<NTHREADS; i++) 
			pthread_join(thr[i],NULL);

//	    	end1 = clock ();
//    		end2 = time (NULL);

  	  	gettimeofday(&tv2, NULL);
 
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
/*
	//Printing the graph
	node *ptr;
	for(int i=1; i<=numver; i++) 
	{
		cout<<ver[i].Vid<<": ";
		ptr=ver[i].next;
		while(ptr!=NULL) 
		{
			cout<<ptr->id<<",";
			cout<<ptr->color<<" ";
			ptr=ptr->next;
		}
		cout<<"Partition: "<<ver[i].Partition_id<<" ";
		cout<<"Boundary vertex: "<<ver[i].boundary_vertex;
		cout<<" Color: "<<ver[i].color;
		cout<<endl;
	}
	cout<<"Max Degree: "<<max_deg<<endl;
*/
	outfile<<"No of Threads : "<<NTHREADS<<endl;

//    	cout << "Duration 1 (clock function): " << fixed << setprecision (DIGITS) << static_cast<double> ((end1 - start1) / CLOCKS_PER_SEC) << " secs." << endl;
//    	cout << "Duration 2 (time function): " << (end2 - start2) << " secs." << endl;
//    	cout << "Duration 3 (gettimeofday() function): " << difference->secs << "." << difference->usecs<<" secs."<<endl;

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
