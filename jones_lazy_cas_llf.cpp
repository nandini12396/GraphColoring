//Barrier - Varying no of threads (Graph is fixed)
//Input graph must be from 1, ..., n

/*
lazy approach doesnt make sense here because there is no remove operation, so no marked field needed
insert in last - so lock the last set_node and added set_node

non blocking will be simple, no remove, hence no marked field again
*/

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

fstream myfile("b_jp.txt", std::ios_base::in);
ofstream outfile("new_results.txt",std::ios_base::app);

int numver, max_deg=0, p;

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
	int ro;
	int n_wait;
	int degree;
//	bool colored;
	list<int> send_queue;
};

Vertex* ver;

typedef struct set_node
{
	int id;
	int color;	
	bool looked;
//	int marked;	//for remove operation
	#ifdef USE_LAZY
	pthread_mutex_t lock;
	#endif
	struct set_node *next;
}set_node;

set_node **head;

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
		ver[i].ro = -1;
		ver[i].n_wait = 0;
		ver[i].degree = 0;		
//		ver[i].colored = false;
//		ver[i].send_queue = list<int>();	
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
		ver[i].degree = d;
		if(d > max_deg)
			max_deg = d;
	//	ver[i].ro = log2(d);	//LLF

//		ptr = ver[i].next;
//		while(ptr!=NULL)
//		{
//			if(ver[ptr->id].ro == ver[i].ro)
//			{
				ver[i].ro = numver - i;
//				break;
//			}
//			ptr = ptr->next;
///		}
//		int k = 1;
//		while(k <=numver)
//		{
//			if(ver[k].ro == ver[i].ro && k!=i)
//			{
		//		ver[i].ro = rand()%numver;
//				k = 1;
//				continue;
//			}
//			else
//				k++;			
//		}
	}

	outfile<<"Max Degree: "<< max_deg << endl;
}

void* add_set(int k, int v, int c)	//it should return whether it worked, add at last
{
	set_node *new_set_node = new set_node();
	new_set_node->id = v;
	new_set_node->color = c;
	new_set_node->looked = false;
	new_set_node->next = NULL;

	#ifdef USE_LAZY
	pthread_mutex_init(&new_set_node->lock, NULL);
	#endif

	#ifdef USE_LAZY

	while(true)
	{
		set_node *curr = head[k];
		while(curr->next !=NULL)
			curr = curr->next;
	
		pthread_mutex_lock(&new_set_node->lock);
		pthread_mutex_lock(&curr->lock);
		if(curr->next == NULL)		//if not NULL then it means that someone else has inserted in the mean time u grabbed lock, so try again
		{
			curr->next = new_set_node;
		//	cout << "Inserted " << new_set_node->id << endl;
	
			pthread_mutex_unlock(&new_set_node->lock);
			pthread_mutex_unlock(&curr->lock);	
			return curr;	
		}
		pthread_mutex_unlock(&curr->lock);	
		pthread_mutex_unlock(&new_set_node->lock);
	}
//	return curr;	

	#elif USE_CAS

	while(true)
	{
		set_node *curr = head[k];
		while(curr->next !=NULL)
			curr = curr->next;

		if (__sync_bool_compare_and_swap(&curr->next, NULL, (set_node*)new_set_node))
		{
		//	cout << "Inserted " << new_set_node->id << endl;
			return curr;
		}	
	}
	#endif
}		

set_node* iterate(set_node *index)	//pass where the current search was
{
	set_node *curr = index;	
	curr = curr->next;
	return curr;
}

void* color(void* t) 
{
	long tid=(long)t;

	list<Vertex*> ulist;     //all vertices of this thread
//	list<Vertex*> recolor;

	int k;
	int avclr[max_deg + 1];

	for(int j=0; j < max_deg+1; j++)  //initialise available array
		avclr[j]=j;

	long partition_id;

	for(int i=1;i<=numver;i++)
	{
		partition_id = (long)ver[i].Partition_id;
		if(partition_id == tid)
		{
//			cout << i;
			ulist.push_back(&ver[i]);
		}
	}

//	cout << tid << endl;
	int temp1,ret,iteration = 0;
	list<int> t1;
	int i=0,size;
//	size = ulist.size();
	node *tp,*tp1;
	int w;

	if(ulist.size() == 0)		
	{
		cout<<"ERROR: "<< tid << " has ulist initially empty!" <<endl;
		exit(0);
	}

	list<int> color_queue;
	int n_colored, n_boundary = 0;

	color_queue.clear();

	node *ptr;
/*        for(int i=1; i<=numver; i++)
        {
                cout<<ver[i].Vid<<": ";
                ptr=ver[i].next;
                while(ptr!=NULL)
                {
                        cout<<ptr->id<<",";
//                      cout<<ptr->color<<" ";
                        ptr=ptr->next;
                }
                cout<<"Partition: "<<ver[i].Partition_id<<" ";
                cout<<"Boundary vertex: "<<ver[i].boundary_vertex << " ";
//                cout<<"Ro: "<<ver[i].ro;
                cout<<" Color: "<<ver[i].color;
                cout<<endl;
        }
*/
	for(list<Vertex*>::iterator iter=ulist.begin(); iter!=ulist.end(); iter++) 
	{
		if((*iter)->color != -1)
			cout << "yes, gone mad" << endl;
		if((*iter)->boundary_vertex == true)	// boundary vertices
		{
//			cout << (*iter)->Vid << " ";
			(*iter)->n_wait = 0;
			(*iter)->send_queue.clear();

			tp = (*iter)->next;
			while(tp!=NULL)
			{
				if(ver[tp->id].boundary_vertex == true)
				{
					w = tp->id;
					if(ver[w].ro > (*iter)->ro)
						(*iter)->n_wait++;
					else
						(*iter)->send_queue.push_back(w);						
				}
				tp = tp->next;
			}
			if((*iter)->n_wait == 0)
				color_queue.push_back((*iter)->Vid);			

			n_boundary++;

//			if((*iter)->degree != (*iter)->n_wait + (*iter)->send_queue.size())	//because of internal vertices this is false claim
//				cout << "gone mad" << endl;
		}
	}

	for(list<int>::iterator it=color_queue.begin(); it!=color_queue.end(); it++) 	//seq-color
	{
		ver[*it].color = 0;
	}
	
	n_colored = color_queue.size();
	
	list<int> temp_list;

	for(list<int>::iterator it=color_queue.begin(); it!=color_queue.end(); it++) 
	{
//		cout << " n wait of color queue " << ver[*it].n_wait << endl;
		temp_list.clear();

		for(list<int>::iterator itr=ver[*it].send_queue.begin(); itr!=ver[*it].send_queue.end(); itr++) 
		{
			temp_list.push_back(ver[*itr].Partition_id);
		}
		temp_list.sort();
		temp_list.unique();

		for(list<int>::iterator itr=temp_list.begin(); itr!=temp_list.end(); itr++)
		{
			k = *itr;
			add_set(k, (*it), ver[*it].color);
		}
	}	
	
	color_queue.clear();

	i = tid;
	set_node *curr = head[i];	//skip sentinel node
	set_node *temp;

	while(n_colored < n_boundary)          //make sure initially ulist has at least 1 vertex
	{
		while(curr->next == NULL)
		{
		}

		temp = iterate(curr);
		if(temp != NULL)
			curr = temp;

	//	k = ver[curr->id].Partition_id;
		
		if(curr->looked == false && curr->id != -1)	//set contains vertex which is colored
		{
			tp = ver[curr->id].next;
			while(tp != NULL)
			{
				if(ver[tp->id].Partition_id == tid && ver[tp->id].ro < ver[curr->id].ro && ver[tp->id].boundary_vertex == true && ver[tp->id].color == -1)	//neighbours of same partition
				{
					w = tp->id;
					tp1 = ver[w].next;		//mark curr's color as forbidden for its neighbour of Ti's partition
					while(tp1 != NULL)
					{
						if(tp1->id == curr->id)
						{
							tp1->color = curr->color;
							break;
						}
						tp1 = tp1->next;
					}

 //                              if(ver[tp->id].Partition_id == tid && ver[tp->id].ro < ver[curr->id].ro)        //neighbours of same partition
 //                               {
					ver[w].n_wait--;
					if(ver[w].n_wait == 0)
						color_queue.push_back(w);
				}
				tp = tp->next;
			}
			curr->looked = true;
		}

		for(list<int>::iterator it=color_queue.begin(); it!=color_queue.end(); it++) 	//seq-color
		{
			t1.clear();
			tp = ver[*it].next;
			while(tp != NULL)
			{
				if(tp->color != -1)
					t1.push_back(tp->color);
//					t1.push_back(ver[tp->id].color);
				tp = tp->next;		
			}			
			k=0;
			list<int>::iterator ti = t1.begin();
			while(ti != t1.end())
			{
				if(*ti == k)
				{
					k++;
					ti = t1.begin();
					continue;
				}
				ti++;
			}
			ver[*it].color = k;
//		        ver[*it].colored = true;

			tp = ver[*it].next;
			while(tp != NULL)
			{
				if(tp->color == k)
					cout << "Something wrong! " << tp->color << endl;
				tp = tp->next;
			}

			tp = ver[*it].next;
	                while(tp != NULL)
        	        {
//				if(ver[tp->id].colored == false)
//				{
				if(ver[tp->id].Partition_id == tid)
				{				
	                        	tp1 = ver[tp->id].next;
//				cout << tp1->id << " " << ver[*it].Vid << endl;
		                        while(tp1 != NULL)
        		                {
						if(tp1->id == *it)
						{
							tp1->color = k;
							break;
						}
					        tp1 = tp1->next;
					}
				}
//				}
//                	        if(tp1->id == *it)
  //                      	        tp1->color = k;
	                        tp = tp->next;
        	        }
		}

		n_colored = n_colored + color_queue.size();

        for(list<int>::iterator it=color_queue.begin(); it!=color_queue.end(); it++)
        {
//              cout << " n wait of color queue " << ver[*it].n_wait << endl;
                temp_list.clear();

                for(list<int>::iterator itr=ver[*it].send_queue.begin(); itr!=ver[*it].send_queue.end(); itr++)
                {
                        temp_list.push_back(ver[*itr].Partition_id);
                }
		temp_list.sort();        
	        temp_list.unique();

                for(list<int>::iterator itr=temp_list.begin(); itr!=temp_list.end(); itr++)
                {
                        k = *itr;
                        add_set(k, (*it), ver[*it].color);
                }
         }


/*
		for(list<int>::iterator it=color_queue.begin(); it!=color_queue.end(); it++) 	
		{
			for(list<int>::iterator itr=ver[*it].send_queue.begin(); itr!=ver[*it].send_queue.end(); itr++) 
			{
				k = ver[*itr].Partition_id;
				add_set(k, (*it), ver[*it].color);
			}	
		}
*/
		color_queue.clear();

//		cout << n_boundary - n_colored << endl;
	}

	for(list<Vertex*>::iterator iter=ulist.begin(); iter!=ulist.end(); iter++) 
	{
		if((*iter)->boundary_vertex != true)	// internal vertices
		{
			t1.clear();
			tp = (*iter)->next;
			while(tp != NULL)
			{
				if(ver[tp->id].color != -1)
					t1.push_back(ver[tp->id].color);
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

                        tp = (*iter)->next;
                        while(tp != NULL)
                        {
                                if(tp->color == (*iter)->color)
                                        cout << "Something mad! " << tp->color << endl;
                                tp = tp->next;
                        }


		}
//		else
//			cout << " n_wait " << (*iter)->n_wait << endl;
	}
}

int main(int argc,char *argv[]) 
{
	myfile>>numver;
	cout<<"The no of vertices: "<<numver<<endl;
	outfile<<"orkut DATA"<<endl;
	outfile<<"USING JONES"<<endl;
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

	ver = new Vertex[numver+1];

	input();


	int i,j,no_colors=0,id=0;

	head = (set_node**) malloc(sizeof(set_node*) * NTHREADS);

	for(i=0;i<NTHREADS;i++)
	{
		head[i] = (set_node*)malloc(sizeof(set_node));		//set Si is represented by head[i]	
//		head[i] = new set_node;
		head[i]->id = -1;	//sentinel set_node
		head[i]->color = -1;
		head[i]->next = NULL;
	
		#ifdef USE_LAZY
		pthread_mutex_init(&head[i]->lock, NULL);
		#endif	
	}

	int temp,dig;
	node* ptr;

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
//	{
		cnt = 0;

		for(i=1;i<=numver;i++)
		{
			ver[i].color = -1;
			ver[i].Partition_id = -1;
			ver[i].boundary_vertex = false;
			ptr = ver[i].next;
			while(ptr != NULL)
			{
				ptr->color = -1;
				ptr = ptr->next;
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

  /*      for(int i=1; i<=numver; i++)
        {
                cout<<ver[i].Vid<<": ";
                ptr=ver[i].next;
                while(ptr!=NULL)
                {
                        cout<<ptr->id<<",";
//                      cout<<ptr->color<<" ";
                        ptr=ptr->next;
                }
                cout<<"Partition: "<<ver[i].Partition_id<<" ";
                cout<<"Boundary vertex: "<<ver[i].boundary_vertex << " ";
	                cout<<"Ro: "<<ver[i].ro;
                cout<<" Color: "<<ver[i].color;
                cout<<endl;
        }
*/
	int i=1,var =0,j;
	while(i <=  numver)
	{
		for(j=0; j<p && i<=numver; j++)
		{
			if(var == NTHREADS - 1)
			{
				while(i<=numver)
				{
					ver[i].Partition_id = var;
					i++;
//					break;
				}
			}
			else
			{
				ver[i].Partition_id = var;
				i++;
			}
		}
		var++;

	}

	for(i=1;i<=numver;i++)
		if(ver[i].Partition_id == -1)
			ver[i].Partition_id = NTHREADS-1;
/*	
	for(i=1;i<=numver;i++)
               if(ver[i].Partition_id != ver[i+1].Partition_id && i != numver)
               {
        	       cout << i << " " << ver[i].Partition_id << endl;
                       cout << i+1 << " "<<ver[i+1].Partition_id << endl;
               }
*/
// node *ptr;
/*        for(int i=1; i<=numver; i++)
        {
                cout<<ver[i].Vid<<": ";
                ptr=ver[i].next;
                while(ptr!=NULL)
                {
                        cout<<ptr->id<<",";
//                      cout<<ptr->color<<" ";
                        ptr=ptr->next;
                }
                cout<<"Partition: "<<ver[i].Partition_id<<" ";
                cout<<"Boundary vertex: "<<ver[i].boundary_vertex << " ";
//                cout<<"Ro: "<<ver[i].ro;
                cout<<" Color: "<<ver[i].color;
                cout<<endl;
       }
*/
			for(i=1; i<=numver; i++) 
			{
				ptr=ver[i].next;
				while(ptr!=NULL) 
				{
					temp=ptr->id;
					if(ver[temp].Partition_id!=ver[i].Partition_id) 
					{
		     				ver[i].boundary_vertex=true;
						break;
					}
					ptr=ptr->next;		
				}
			}
	//		cout << "boundary done";
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
/*       node *ptr;
        for(int i=1; i<=numver; i++)
        {
                cout<<ver[i].Vid<<": ";
                ptr=ver[i].next;
                while(ptr!=NULL)
                {
                        cout<<ptr->id<<",";
//                      cout<<ptr->color<<" ";
                        ptr=ptr->next;
                }
                cout<<"Partition: "<<ver[i].Partition_id<<" ";
                cout<<"Boundary vertex: "<<ver[i].boundary_vertex << " ";
                cout<<"Ro: "<<ver[i].ro;
                cout<<" Color: "<<ver[i].color;
                cout<<endl;
       }

*/

		gettimeofday(&tv1,NULL);
	
		cout<<"calling func"<<endl;
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
//	}

	double average = accumulate(vec.begin(),vec.end(),0.0) / vec.size();

//	node *ptr;
/*	for(int i=1; i<=numver; i++) 
	{
		cout<<ver[i].Vid<<": ";
		ptr=ver[i].next;
		while(ptr!=NULL) 
		{
			cout<<ptr->id<<",";
//			cout<<ptr->color<<" ";
			ptr=ptr->next;
		}
		cout<<"Partition: "<<ver[i].Partition_id<<" ";
		cout<<"Boundary vertex: "<<ver[i].boundary_vertex << " ";
		cout<<"Ro: "<<ver[i].ro;
		cout<<" Color: "<<ver[i].color;
		cout<<endl;
	}

*/

	//Printing the graph
//	node *ptr;
/*	for(int i=1; i<=numver; i++) 
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
	int fal = 0;
  /*     for(i=1;i<=numver;i++)
        {
                if(ver[i].color == -1)
                       cout << "ERROR: node with id " << i << " not getting colored at all";
                ptr = ver[i].next;
                while(ptr!=NULL)
                {
                        if(ver[i].color == ver[ptr->id].color)
//				break;	
		cout << i << " of partition " << ver[i].Partition_id <<  " and " << ptr->id << " of partition " << ver[ptr->id].Partition_id << " have same color." << endl; // cout << "Graph is not properly colored.";
                        ptr = ptr->next;
                }
		if(ptr!=NULL)
			fal++;
        }
	cout << "No of vertices not colored properly: " << fal;

	for(int i=1; i<=numver; i++) 
	{
		cout<<ver[i].Vid<<": ";
		ptr=ver[i].next;
		while(ptr!=NULL) 
		{
			cout<<ptr->id<<",";
//			cout<<ptr->color<<" ";
			ptr=ptr->next;
		}
		cout<<"Partition: "<<ver[i].Partition_id<<" ";
		cout<<"Boundary vertex: "<<ver[i].boundary_vertex << " ";
		cout<<"Ro: "<<ver[i].ro;
		cout<<" Color: "<<ver[i].color;
		cout<<endl;
	}
*/
      for(i=1;i<=numver;i++)
        {
                if(ver[i].color == -1)
                       cout << "ERROR: node with id " << i << " not getting colored at all";
                ptr = ver[i].next;
                while(ptr!=NULL)
                {
                        if(ver[i].color == ver[ptr->id].color)
                        {
				fal++;
				 cout << i << " of partition " << ver[i].Partition_id << " " << ver[i].boundary_vertex <<  " and " << ptr->id << " of partition " << ver[ptr->id].Partition_id << " " << ver[ptr->id].boundary_vertex << " have same color " << ver[i].color << endl; // cout << "Graph is not properly colored.";
			}
                           ptr = ptr->next;
                 }
 //                if(ptr!=NULL)
   //                     fal++;
         }
//
	cout << "No of false = " << fal << endl;
	outfile << "No of false = " << fal << endl;
//	free(ver);
//	free(difference);
//	free(thr);
	
	return 0;
}
