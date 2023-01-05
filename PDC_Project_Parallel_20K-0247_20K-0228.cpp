// CPP program to implement sequence alignment
// problem.
#include <bits/stdc++.h>
#include <omp.h>
#include <vector>

int total = 0;

using namespace std;

int Find_Minimum_Of_All_Three_Options(int a, int b, int c) 
{
	if (a <= b && a <= c) 
    {
		return a;
	}
	else if (b <= a && b <= c) 
    {
		return b;
	}
	else 
    {
		return c;
	}
}

//Task3 function for 1st section (description given later). Here the only concern while applying parallel for clause is of updating shared variable i and using
//it as index to fill up the xans vector. Hence it was required to either lock the variable or the code snippet or apply reduction. Locking the code snippet is 
//not required here so no use of critical as it is less efficient also. Also we can apply atomic, but what if hardware instructions are not available at the time of code execution. Therefore, we thought to use
//reduction.

void Task3(int xpos, int i, vector<int> &xans, vector<char> &x, int numthreads){
	#pragma omp parallel for reduction(-:i) num_threads(numthreads)
		for (int k = xpos; k > 0; k--)
		{
			if (i > 0)
			{
				i-=1;
				xans[k] = (int)x[i];
			} 
			else xans[k] = (int)'_';
		}
}

//Task4 function for 2nd section (description given later). Here the only concern while applying parallel for clause is of updating shared variable j and using
//it as index to fill up the yans vector. Hence it was required to either lock the variable or the code snippet or apply reduction. Locking the code snippet is 
//not required here so no use of critical as it is less efficient also. Also we can apply atomic, but what if hardware instructions are not available at the time of code execution. Therefore, we thought to use
//reduction.

void Task4(int ypos, int j, vector<int> &yans, vector<char> &y, int numthreads){
	#pragma omp parallel for reduction(-:j) num_threads(numthreads)
		for (int k = ypos; k > 0; k--)
		{
			if (j > 0)
			{
				j-=1;
				yans[k] = (int)y[j];
			}
			else yans[k] = (int)'_';
		}
}

// function to find out the minimum penalty
void getMinimumPenalty(vector<char> &x, vector<char> &y, int pxy, int pgap, int numthreads)
{
	
	int m = x.size(); // length of gene1
	int n = y.size(); // length of gene2
	
	// table for storing optimal substructure answers

	vector<vector<int> > dp(m+1);
	
	#pragma omp parallel for num_threads(numthreads)
		for (int i = 0; i < m+1; i++)
		{
			dp[i] = vector<int>(n+1);
		}

	//Simply parallelism applied on for loop without any concerns as the shared variable initialization is independent of each other so we can distribute it amongst different number of threads.

	#pragma omp parallel for num_threads(numthreads)
		for (int i = 0; i < n+1; i++)
		{
			dp[0][i] = i * pgap;
		}

	//Simply parallelism applied on for loop without any concerns as the shared variable initialization is independent of each other so we can distribute it amongst different number of threads.
	
	#pragma omp parallel for num_threads(numthreads)
		for(int i= 0; i < m+1; i++)
		{
			dp[i][0] = i * pgap;
		}

	// calculating the minimum penalty
	int End, Start;

	if (m == n)
	{
		Start = 1;
		End = n;
	}
	else if (m < n)
	{
		End = m; 
		Start = 0;
	}
	else
	{
		End = n;
		Start = m - n + 1;
	}

	//The code (of filling out dp matrix) is modified to achieve parallelism for the very important part of code. The original block of code was not supporting the 
	//parallelism. The parallelism for the program except that part of code was not worth.

	//To calculate the matrix dp as per the rules of dynamic programming, it was observed that at each iteration of nested for loop, the previous elements are re-
	//quired. Hence it could be dangerous to use parallel for clause for outer loop as well as inner loop. Therefore, we decided to go for inner loop parallelism.

	#pragma omp parallel num_threads(numthreads)
	for (int k = 2; k <= End; k++)
	{
		#pragma omp for
		for (int j = k - 1; j > 0; j--)
		{
			
			int i = k - j;
			if (x[i - 1] == y[j - 1])
			{
				dp[i][j] = dp[i - 1][j - 1];
			}
			else
			{
				dp[i][j] = Find_Minimum_Of_All_Three_Options(dp[i - 1][j - 1] + pxy,
					dp[i - 1][j] + pgap,
					dp[i][j - 1] + pgap);
			}
		}
	}

	//To calculate the matrix dp as per the rules of dynamic programming, it was observed that at each iteration of nested for loop, the previous elements are re-
	//quired. Hence it could be dangerous to use parallel for clause for outer loop as well as inner loop. Therefore, we decided to go for inner loop parallelism.

	if (m < n)
	{
		#pragma omp parallel num_threads(numthreads) 
		for (int k = m + 1; k <= n; k++)
		{
			#pragma omp for 
			for (int i = 1; i <= m; i++)
			{
				int j = k - i;
				if (x[i - 1] == y[j - 1])
				{
					dp[i][j] = dp[i - 1][j - 1];
				}
				else
				{
					dp[i][j] = Find_Minimum_Of_All_Three_Options(dp[i - 1][j - 1] + pxy,
					dp[i - 1][j] + pgap,
					dp[i][j - 1] + pgap);
				}
			}
		}
	}

	//To calculate the matrix dp as per the rules of dynamic programming, it was observed that at each iteration of nested for loop, the previous elements are re-
	//quired. Hence it could be dangerous to use parallel for clause for outer loop as well as inner loop. Therefore, we decided to go for inner loop parallelism.

	else if (m > n)
	{
		#pragma omp parallel num_threads(numthreads)
		for (int k = n + 1; k <= m; k++)
		{
			#pragma omp for 
			for (int j = n; j > 0; j--)
			{

				int i = k - j;
				if (x[i - 1] == y[j - 1])
				{
					dp[i][j] = dp[i - 1][j - 1];
				}
				else
				{
					dp[i][j] = Find_Minimum_Of_All_Three_Options(dp[i - 1][j - 1] + pxy,
					dp[i - 1][j] + pgap,
					dp[i][j - 1] + pgap);
				}
			}
		}
	}

	if (Start == 0)
	{
		Start++;
	}

	//To calculate the matrix dp as per the rules of dynamic programming, it was observed that at each iteration of nested for loop, the previous elements are re-
	//quired. Hence it could be dangerous to use parallel for clause for outer loop as well as inner loop. Therefore, we decided to go for inner loop parallelism.


	#pragma omp parallel num_threads(numthreads)
	for (int k = Start; k <= m; k++)
	{
		#pragma omp for
		for (int i = k; i <= m; i++)
		{
			int j = k + n - i;
			if (i != 0 || j != 0)
			{
				if (x[i - 1] == y[j - 1])
				{
					dp[i][j] = dp[i - 1][j - 1];
				}
				else
				{
					dp[i][j] = Find_Minimum_Of_All_Three_Options(dp[i - 1][j - 1] + pxy,
					dp[i - 1][j] + pgap,
					dp[i][j - 1] + pgap);
				}
			}
		}

	}

	// Reconstructing the solution
	int l = n + m; // maximum possible length
	
	int i = m, j = n;
	
	int xpos = l;
	int ypos = l;

	// Final answers for the respective strings
	vector<int> xans(l+1,0);
	vector<int> yans(l+1,0);

	//For the while loop, it can easily be converted in for loop, we tried to do so in order to apply parallel for clause. Since each and every variable is being
	//updated at each iteration, it would require to block whole section of code (if else if clauses) OR apply multiple critical sections to synchronize the upda-
	//tion of each variable. Hence being time consuming and increasing parallel overhead, it is better to sequentially execute this while loop.

	while(!(i == 0 || j == 0))
	{
		if (x[i - 1] == y[j - 1])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j - 1] + pxy == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)'_';
			i--;
		}
		else if (dp[i][j - 1] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)'_';
			yans[ypos--] = (int)y[j - 1];
			j--;
		}
	}
	
	//We were trying to nest for loop within sections. To complete the initialization of arrays xans and yans, being independent, we created two sections for them.
	//Later we come to know that we can't nest for loop within sections as for loop require all threads to participate in its iterations whereas the threads are
	//divided between sections. Hence we thought to make 2 functions, each called by respective section. Then we used parallel for clause within those functions.
	
	#pragma omp parallel sections num_threads(2)
	{
		#pragma omp section
		{
			Task3(xpos,i,xans,x, numthreads);
		}
		#pragma omp section
		{
			Task4(ypos,j,yans,y, numthreads);
		}
	}

	// Since we have assumed the answer to be n+m long,
	// we need to remove the extra gaps in the starting
	// id represents the index from which the arrays
	// xans, yans are useful
	int id = 1;
	
	//We modified the working of for loop. The original code had for loop executing in reverse. We tried to execute it in forward direction with parallel for
	//clause so that we can store the last index in id variable where the condition met. But still the the iteration of k that we require to get stored in id
	//is possible to be executed before. CORRECT OUTPUT IS BETTER THAN FALSE PARALLELISM.
	
	//#pragma omp parallel for num_threads(numthreads)
		for (int k = 1; k <= l; k++)
		{
			if ((char)yans[k] == '_' && (char)xans[k] == '_')
			{
				id = k + 1;
			}
		}

	// Printing the final answer
	total = total + dp[m][n];

	//To store the updated aligned strings, we made sections. As we know that storing two different strings in two different files are completely independent tasks.
	//Less benefit is better than 0. Hence we can apply sections for them.


	#pragma omp parallel sections num_threads(2)
	{
		#pragma omp section
		{
			fstream file;
			file.open("gene1.txt",ios_base::app);
			for (int k = id; k <= l; k++)
			{
				if(xans[k]>0)
				{
					file.put((char)xans[k]);	
				}
			} 
			 
			file.close();
		}
		#pragma omp section
		{
			fstream file;
			file.open("gene2.txt",ios_base::app);
			for (int k = id; k <= l; k++)
			{
				if(yans[k]>0)
				{
					file.put((char)yans[k]);	
				}
			}
			 
			file.close();
		}
	}

	return;
}

// Driver code
int main(int argc, char *argv[]){
    

    // Required code for which execution time needs to be computed
	int numthreads = atoi(argv[1]);

	// input strings
	string gene1, gene2;
	int misMatchPenalty, gapPenalty;

	cout<<"Enter string 1: ";
	cin>>gene1;
	cout<<"Enter string 2: ";
	cin>>gene2;
	cout<<"Enter the mis matching penalty: ";
	cin>>misMatchPenalty;
	cout<<"Enter the gap penalty: ";
	cin>>gapPenalty;
	double itime, ftime, exec_time;
	int j,k;
	itime = omp_get_wtime();
	for (int i = 0; i < 500; i++) //execute it 100 times to concatenate for 1 million, 300 time for 3 million, and 500 times for 5 million
	{	
		vector<char> genee1(10000,0);
		vector<char> genee2(10000,0);
	
		j=0;
		k=0;
		for (int i = 0; i < 10000; i++)
		{
			genee1.push_back(gene1[j]);
			genee2.push_back(gene2[k]);
			if (j == int(gene1.length()-1))
			{
				j = 0;
			}
			if (k == int(gene2.length()-1))
			{
				k = 0;
			}
			j++;
			k++;
		}

		// calling the function to calculate the result
		getMinimumPenalty(genee1, genee2, misMatchPenalty, gapPenalty, numthreads);
	}
	ftime = omp_get_wtime();
	exec_time = ftime - itime;
	cout << "The aligned genes are stored in file ""gene1.txt"" and ""gene2.txt"":\n";
	cout << "Minimum Penalty in aligning the genes = ";
	cout << total << "\n";
	cout << "\nTime taken is: "<< exec_time<<"\n";
	fstream file;
	file.open("ExecutionDetailsLogFile_Parallel_20K-0247_20K-0228.txt",ios_base::app);
	file<<"The time takein by parallel code for two strings is: "<<exec_time<<" on "<<numthreads<<" threads.\n";
	file.close();
	return 0;
}
