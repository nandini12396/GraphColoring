// BTO/SGT Protocol of STMs - Varying no of threads (Graph is fixed)
// NOTE: Input graph should have vertices from 1...n (should not have 0 as a vertex)
// g++ -std=c++11 -DUSE_BTO tm_graph.cpp -lpthread -ltbb
// g++ -std=c++11 -DUSE_SGT tm_graph.cpp -lpthread -ltbb

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

#ifdef USE_BTO
#include "BTO.cpp"
#elif USE_SGT
#include "SGT.cpp"
#elif USE_MVTO
#include "MVTO.cpp"
#endif 

using namespace std;

#define DIGITS 4

fstream myfile("input.txt", std::ios_base::in);
ofstream outfile("paper.txt",std::ios_base::app);

int numver, max_deg=0, p;

#ifdef USE_BTO
STM *lib = new BTO;
#elif USE_SGT
STM *lib = new SGT;
#elif USE_MVTO
STM *lib = new MVTO(1);
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
};

Vertex* ver;

int aborts = 0;

pthread_mutex_t lock1;  //for updating aborts

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
		ver[i].next = NULL;
		ver[i].Vid = i;
		ver[i].boundary_vertex = false;
		ver[i].color = -1;
		ver[i].Partition_id = -1;
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
	outfile<<"No of undirected edges: "<< edge << endl;
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

	list<Vertex*> ulist; //all internal vertices of this thread's partition
	list<Vertex*> blist; //boundary vertices of this thread's partition

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
//		temp.sort();
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
//	for(list<Vertex*>::iterator iter=ulist.begin(); iter!=ulist.end(); iter++) 
//	{
//		cout<<(*iter)->Vid<<" "<<(*iter)->color<<endl;
//	}

	cout<<tid<<"'s Internal Vertices Colored"<<endl;

    	common_tOB* tob1=new common_tOB;
    	tob1->value=operator new(sizeof(Vertex));   /*for int*/
    	tob1->size=sizeof(Vertex);

	int flag; //return value
	int num_abort = 0;
	long long err;
	Vertex v;

	//write the colors of internal vertices to the shared memory graph

	for(list<Vertex*>::iterator iter=ulist.begin(); iter!=ulist.end(); iter++) 
	{
l2:		trans_state *T = lib->begin(); 
		tob1->ID = (*iter)->Vid;
		
		v.color = (*iter)->color;
		v.boundary_vertex = (*iter)->boundary_vertex;
		v.Vid = (*iter)->Vid;
		v.Partition_id = (*iter)->Partition_id;
		v.next = (*iter)->next;

		*((Vertex*)tob1->value) = v;

		flag = 0;

		flag = lib->write(T,tob1);		//write the whole vertex again
		if(flag != 0)
			cout << "write fail" << endl;	

		flag = lib->try_commit(T,err);		//here the no of aborts should be 0, each writing to diff location
		if(flag != 0)
		{
			num_abort++;
			goto l2;
		}
	}

	for(list<Vertex*>::iterator iter=blist.begin(); iter!=blist.end(); iter++) 
	{
//		cout<<(*iter)->Vid<<" getting colored"<<endl;
l1:		trans_state *T = lib->begin();         /* coloring of each boundary vertex is handled by a transaction*/
				
		tp=(*iter)->next;
		t1.clear();

		//Check the colors of adjacent vertices and assign a least possible color different from the adjacent vertices.
		while(tp!=NULL) 
		{
			tob1->ID = tp->id;
			flag = lib->read(T,tob1);
			if(flag != 0)
			{
				num_abort++;
				goto l1;				
			}
//			cout<<(*(Vertex*)tob1->value).color<<endl;
			if((*(Vertex*)tob1->value).color!=-1)
				t1.push_back((*(Vertex*)tob1->value).color);
			tp=tp->next;
		}
		k=0;
//		t1.sort();
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

		tob1->ID = (*iter)->Vid; 	//coloring the boundary vertex
		
		v.color = k;
		v.boundary_vertex = (*iter)->boundary_vertex;
		v.Vid = (*iter)->Vid;
		v.Partition_id = (*iter)->Partition_id;
		v.next = (*iter)->next;

		*((Vertex*)tob1->value) = v;

		flag = -1;

	//	while(flag != 0)
//		{
			flag = lib->write(T,tob1);		//write the whole vertex again
			if(flag != 0)
				cout << "write fail" << endl;	
	//			num_abort++;
//		}

		flag = -1;

//		while(flag != 0)
//		{
			flag = lib->try_commit(T,err);		//write the whole vertex again
			if(flag != 0)
			{
				num_abort++;
				goto l1;		
			}
//		}

		(*iter)->color = k;
//		cout<<(*iter)->Vid<<" "<<(*iter)->color<<endl;
	}

//	trans_state *T = lib->begin();
//	for(list<Vertex*>::iterator iter=blist.begin(); iter!=blist.end(); iter++)
	{
//		tob1->ID = (*iter)->Vid;
//		lib->read(T,tob1);
//		cout<<(*(Vertex*)tob1->value).Vid<<" "<<(*(Vertex*)tob1->value).color<<endl;	
	}
//	lib->try_commit(T,err);
	
//	for(list<Vertex*>::iterator iter=blist.begin(); iter!=blist.end(); iter++) 
//	{
//		cout<<(*iter)->Vid<<" "<<(*iter)->color<<endl;
//	}

	cout<<tid<<"'s aborts = "<<num_abort<<endl;
	outfile<<tid<<"'s aborts = "<<num_abort<<endl;

	pthread_mutex_lock(&lock1);
	aborts = aborts + num_abort;
	pthread_mutex_unlock(&lock1);

	cout<<tid<<"'s boundary vertices completed.."<<endl;
}

int main(int argc,char *argv[]) 
{
	myfile >> numver;
	cout<<"The no of vertices: "<<numver<<endl;
	
	#ifdef USE_BTO
	outfile<<"Using BTO Protocol:"<<endl;
	#elif USE_SGT
	outfile<<"Using SGT Protocol:"<<endl;
	#elif USE_MVTO
	outfile<<"Using MVTO Protocol:"<<endl;
	#else
	cout<<"Error: Use a protocol (BTO/SGT/MVTO) to compile"<<endl;	
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

	ver=(Vertex*) (malloc((numver+1) * sizeof(Vertex)));

	input();

	int i,j;

	//Creating the graph in shared memory too
	long long err;
	int flag = 0;

	for(i=1;i<=numver;i++)	
		if(lib->create_new(i,sizeof(Vertex)) != SUCCESS)
		{
			cout << "Failed to create id " << i << endl;	
			return 0;
		}

	common_tOB* tob=new common_tOB;
	tob->size=sizeof(Vertex);
	tob->value=operator new(sizeof(Vertex));

	for(i=1;i<=numver;i++)
	{
		trans_state* T=lib->begin();  //Since it is accessed sequentially, it will commit definitely
		tob->ID=ver[i].Vid;
		*((Vertex*)tob->value)=ver[i];
		flag = lib->write(T,tob);
		if(flag == 0)
			lib->try_commit(T,err);		
	}

//	j = lib->try_commit(T,err);
//	if(j==0)
//		cout<<"Shared Memory Graph Created"<<endl;

/*	T = lib->begin();
	for(i=1;i<=numver;i++)
	{
		tob->ID = i;
		lib->read(T,tob);
		cout<<(*(Vertex*)tob->value).Vid<<" created"<<endl;		
	}
*/
	aborts = 0;
	pthread_mutex_init(&lock1,NULL);

	double duration;
 	time_t start1, start2;
    	time_t end1, end2;

	int no_colors=0;

	float x = (float)numver / (float)NTHREADS;
	p=floor(x); //no of vertices in 1 partition

	//Assigning Partition Id's randomly

	srand(unsigned(std::time(0)));
	vector<int> myvector;
	for (int i=1; i<=numver; ++i) 
		myvector.push_back(i);

//	vector<int>::iterator itemp=myvector.begin();
//	while(itemp != myvector.end())
//	{
//		cout<<*itemp<<" ";
//		itemp++;
//	}
//	cout<<"After shuffling";

	random_shuffle ( myvector.begin(), myvector.end() );
	vector<int>::iterator itemp=myvector.begin();
//	itemp = myvector.begin();
//	while(itemp != myvector.end())
//	{
//		cout<<*itemp<<" ";
//		itemp++;
//	}
//	cout<<"No of vertices in 1 thread: "<<p;

			int u=0,z=0;
			i =1;

	  while(i <=  numver)
          {
                  for(j=0; j<p && i<=numver; j++)
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


//	for(i=1;i<=numver;i++)
//			cout<<ver[i].Vid<<" "<<ver[i].Partition_id<<endl;

	//Identifying boundary vertices 

	int temp;
	node* pt;

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
		cout<<"Partition: "<<ver[i].Partition_id<<" ";
		cout<<"Boundary vertex: "<<ver[i].boundary_vertex<<" ";
		cout<<"Color: "<<ver[i].color;
		cout<<endl;
	}
	cout<<"Max Degree: "<<max_deg<<endl;
*/

	pthread_t *thr = new pthread_t [NTHREADS];
    	// Make threads Joinable for sure.
    	pthread_attr_t attr;
    	pthread_attr_init (&attr);
    	pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
 	struct timeval tv1, tv2;
	TIME_DIFF * difference;

	start1 = clock();	
	start2 = time (NULL);
   	gettimeofday(&tv1, NULL);

    	for(int i=0; i<NTHREADS; i++) 
		pthread_create(&thr[i], &attr, color, (void*) i);
	
	for(int i=0; i<NTHREADS; i++) 
		pthread_join(thr[i],NULL);
	
    	end1 = clock ();
    	end2 = time (NULL);
 
   	gettimeofday(&tv2, NULL);
 	difference = my_difftime (&tv1, &tv2);

 	outfile<<"No of Threads : "<<NTHREADS<<endl;
    	
	cout << "Duration 1 (clock function): " << fixed << setprecision (DIGITS) << static_cast<double> ((end1 - start1) / CLOCKS_PER_SEC) << " secs." << endl;
    	cout << "Duration 2 (time function): " << (end2 - start2) << " secs." << endl;
    	cout << "Duration 3 (gettimeofday() function): " << difference->secs << "." << difference->usecs<<" secs."<<endl;

	outfile<<"Time Taken : "<<(end2-start2)<<" secs."<<endl; 
    	outfile << "Duration (gettimeofday() function): " << difference->secs << "." << difference->usecs<<" secs."<<endl;

	no_colors = 0;
	for(i=1;i<=numver;i++)
	{
		if(ver[i].color > no_colors)
			no_colors = ver[i].color;
	}
	no_colors++;

//	for(i=1;i<=numver;i++)
//		cout<<ver[i].Vid<<" "<<ver[i].color<<endl;
	cout<<"No of colors used: " << no_colors<<endl;
	outfile<<"No of colors used: " << no_colors<<endl;
	cout<<"No of aborts: " << aborts<<endl;
	outfile<<"No of aborts: " << aborts<<endl;

	free(ver);
	free(difference);
	free(thr);
	
	return 0;
}
