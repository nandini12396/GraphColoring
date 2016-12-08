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
#include <atomic>

using namespace std;

#define DIGITS 4

fstream myfile("barrier_data.txt", std::ios_base::in);
ofstream outfile("paper.txt",std::ios_base::app);

int numver, max_deg=0, p;

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
	int degree;
	node *next;
};

Vertex *ver;

bool **table;
pthread_mutex_t lock;

atomic<int> *enabled; 
atomic<int> *numReq;
atomic<bool*> head;

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
		ver[i].degree = d;		
		if(d > max_deg)
			max_deg = d;
	}
	outfile<<"Max Degree: "<< max_deg << endl;
}

bool sort_ldf(const Vertex *a, const Vertex *b)
{
	return 	(a->degree) > (b->degree);
}

void request_table(int v, bool** start)
{
	list<int> temp;
	temp.clear();
	node* tp;
	bool flag = false;

	int curr_row = 1, next_row;
	pthread_mutex_lock(&lock);

	*start = head.load();

	while(table[curr_row] != *start)
		curr_row++;	

	tp = ver[v].next;
	temp.push_back(v);
	while(tp != NULL)
	{
		temp.push_back(tp->id);
		tp = tp->next;
	}
//	temp.sort(); no need to sort

	while(true)
	{
		flag = false;
		for(list<int>::iterator t=temp.begin();t!=temp.end();t++)
		{
			if(table[curr_row][*t] == true)
				flag = true;	//go to next row
		}
		if(flag == false)
			break;

		next_row = curr_row + 1;

		curr_row = next_row;
	}

	for(list<int>::iterator t=temp.begin();t!=temp.end();t++)
		table[curr_row][*t] = true;

	numReq[curr_row].fetch_add(1);
		
	pthread_mutex_unlock(&lock);
	
	while(enabled[curr_row].load() != 1)
	{
	}

	*start = table[curr_row];

	return;
}

void release_table(int v, bool** start)
{
	list<int> temp;
	temp.clear();
	node* tp;

	int curr_row = 1, next_row;
	while(table[curr_row] != *start)
		curr_row++;	

	tp = ver[v].next;
	temp.push_back(v);
	while(tp != NULL)
	{
		temp.push_back(tp->id);
		tp = tp->next;
	}
//	temp.sort();

	pthread_mutex_lock(&lock);

	for(list<int>::iterator t=temp.begin();t!=temp.end();t++)
		table[curr_row][*t] = false;

	numReq[curr_row].fetch_sub(1);

	if(numReq[curr_row].load() == 0)
	{
		next_row = curr_row + 1;
		enabled[curr_row].store(0);
		bool *next = table[next_row];
		head.store(next);
		enabled[next_row].store(1); 
	}

	pthread_mutex_unlock(&lock);

	return;
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

//	blist.sort( sort_ldf );
	bool *start;

	for(list<Vertex*>::iterator iter=blist.begin(); iter!=blist.end(); iter++) 
	{
		request_table((*iter)->Vid, &start);

		temp.clear();
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

		release_table((*iter)->Vid, &start);
	}
	cout<<tid<<"'s boundary vertices completed"<<endl;
}

int main(int argc,char *argv[]) 
{
	myfile >> numver;
	cout<<"The no of vertices: "<<numver<<endl;
	outfile << "LIVE JOURNAL DATA" << endl;

	outfile<<"The no of vertices: "<<numver<<endl;
	int NTHREADS, i, j;

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
		
	table = new bool*[numver+1];
	for(int i = 1; i <numver+1; i++)
		table[i] = new bool[numver+1];

	for(i=1;i<=numver;i++)
		for(j=1;j<=numver;j++)
			table[i][j]																																																																																																																																																																																																																																																																																																																																																										 = false;

	enabled = new atomic<int>[numver+1];
	enabled[1].store(1);
	for(i=2;i<=numver;i++)
		enabled[i].store(0);

	numReq = new atomic<int>[numver+1];
	for(i=1;i<=numver;i++)
		numReq[i].store(0);

	head.store(table[1]);

	double duration;
 	time_t start1, start2;
	time_t end1, end2;
   
	int no_colors=0;
	int temp,dig;
	node* pt;
	double b;
	
	cout<<"Graph Generated"<<endl;


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
			p=(int)floor(numver/NTHREADS); //no of vertices in 1 partition
	
			//Assigning Partition Id's randomly
			struct timeval tv3, tv4;
			TIME_DIFF * difference1;

			gettimeofday(&tv3,NULL);
		int u = 0, z = 0;
		i =1;

			 while(i <= numver)
                        {
                                for(j=0;j<p && i<=numver ;j++)
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
					ver[i].Partition_id = NTHREADS - 1;


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
	for(i=1;i<=numver;i++)
	{
		if(ver[i].color == -1)
			cout << "ERROR: node with id " << i << " not getting colored at all";
		pt = ver[i].next;
		while(pt!=NULL)
		{
			if(ver[i].color == ver[pt->id].color)
				cout << "Graph is not properly colored.";
			pt = pt->next;
		}		
	}

	return 0;
}
