// CPP program to implement sequence alignment
// problem.
#include <bits/stdc++.h>
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

void getMinimumPenalty(vector<char> &x, vector<char> &y, int pxy, int pgap)
{
    int i, j; // initialising variables
     
    int m = x.size(); // length of gene1
    int n = y.size(); // length of gene2
     
    // table for storing optimal substructure answers
    vector<vector<int> > dp(m+n+1);
	
    for (int i = 0; i < n+m+1; i++)
    {
		dp[i] = vector<int>(n+m+1);
    }
 
    // initialising the table
    for (i = 0; i <= (n+m); i++)
    {
        dp[i][0] = i * pgap;
        dp[0][i] = i * pgap;
    }   
 
    // calculating the minimum penalty
    for (i = 1; i <= m; i++)
    {
    	
        for (j = 1; j <= n; j++)
        {
        	
            if (x[i - 1] == y[j - 1])
            {
                dp[i][j] = dp[i - 1][j - 1];
            }
            
            else
            {
                dp[i][j] = Find_Minimum_Of_All_Three_Options(dp[i - 1][j - 1] + pxy, dp[i - 1][j] + pgap, dp[i][j - 1] + pgap);
            }
            
        }
        
    }
 
    // Reconstructing the solution
    int l = n + m; // maximum possible length
     
    i = m; j = n;
     
    int xpos = l;
    int ypos = l;
 
    // Final answers for the respective strings
    vector<int> xans(l+1,0);
    vector<int> yans(l+1,0);
     
    while ( !(i == 0 || j == 0))
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
    
    while (xpos > 0)
    {
    	
        if (i > 0)
        {
        	xans[xpos--] = (int)x[--i];
	}
		
        else
        {
        	xans[xpos--] = (int)'_';
	}
		
    }
    
    while (ypos > 0)
    {
    	
        if (j > 0) 
        {
        	yans[ypos--] = (int)y[--j];
	}
		
        else 
        {
        	yans[ypos--] = (int)'_';
	}
		
    }
 
    // Since we have assumed the answer to be n+m long,
    // we need to remove the extra gaps in the starting
    // id represents the index from which the arrays
    // xans, yans are useful
    int id = 1;
    for (i = l; i >= 1; i--)
    {
    	
        if ((char)yans[i] == '_' && (char)xans[i] == '_')
        {
            id = i + 1;
            break;
        }
        
    }
    
    total = total + dp[m][n];
    
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
		
	file.open("gene2.txt",ios_base::app);
	for (int k = id; k <= l; k++)
	{
		
		if(yans[k]>0)
		{
			file.put((char)yans[k]);	
		}
		
	}
			
	file.close();

	return;
}

// Driver code
int main(){
    

    // Required code for which execution time needs to be computed

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
	int j,k;
	time_t start, end;
	time(&start);
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
			if (j == (int)gene1.size()-1)
			{
				j = 0;
			}
			if (k == (int)gene2.size()-1)
			{
				k = 0;
			}
			j++;
			k++;
		}

		// calling the function to calculate the result
		getMinimumPenalty(genee1, genee2, misMatchPenalty, gapPenalty);
	}
	time(&end);
	double time_taken = double(end - start);
	cout << "The aligned genes are stored in file ""gene1.txt"" and ""gene2.txt"":\n";
	cout << "Minimum Penalty in aligning the genes = ";
	cout << total << "\n";
	cout << "The time taken is: " << time_taken  << " seconds\n";
	fstream file;
	file.open("ExecutionDetailsLogFile_Serial_20K-0247_20K-0228.txt",ios_base::app);
	file<<"The time takein by serial code for two strings is: "<<time_taken<<"\n";
	file.close();
	return 0;
}
