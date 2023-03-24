// Jason Clark 
// CS476 Project 2: Phylogenetic Trees

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <sys/stat.h>
#include <iterator>
#include <map>

using namespace std;

struct fileContents;
struct coordinates;

fileContents readFile(string inFile);
void printVector(vector<auto> data);
void printMatrix(vector<vector<auto>> data);
void upgma(fileContents data);
string neighborJoin(fileContents data);
coordinates searchMin(vector<vector<double>> matrix);
vector<double> setPrimeValues(vector<vector<double>> matrix);
vector<vector<double>> setPrimeMatrix(vector<vector<double>> matrix, vector<double> primes);
vector<vector<double>> reduceDistanceMatrix(coordinates coords, vector<vector<double>> matrix, string algo);
vector<pair<string,string>> setClusters(vector<string> data);
vector<pair<string,string>> reduceClusters(coordinates coords, vector<pair<string,string>> data);
void printPrimeValues(vector<pair<string,string>> clusters, vector<double> primes);
void printPrimeMatrix(vector<pair<string,string>> clusters, vector<vector<double>> primeMatrix);

struct fileContents
{
    int size;
    vector<string> nodes;
    vector<vector<double>> matrix;
};

struct coordinates
{
    int x;
    int y;
};

int main()
{
    bool validFile = false;
    string userInput;
    fileContents data;
    data.size = -1;

    cout << "Enter filename: \n>" ; 
    cin >> userInput;
    data = readFile(userInput);
    if(data.size == -1)
        return 0;

    cout << "Enter algorithm to run ('u' = UPGMA, 'n' = Neighbor Join):\n>";
    cin >> userInput;
    if(userInput == "u")
        upgma(data);
    else if(userInput == "n")
        neighborJoin(data);
    else
        cout << "Invalid user input!\n"; 

    cout << endl << "Program has finished execution..." << endl; 
    return 0;
}

void upgma(fileContents data)
{
    double distance;
    vector<pair<string,string>> clusters = setClusters(data.nodes);
    vector<vector<double>> distanceMatrix = data.matrix;
    vector<vector<double>> copy;
    coordinates minCoords;
    printMatrix(distanceMatrix);
    while(clusters.size() > 1)
    {
        //Print Clusters
        for(int i = 0; i < clusters.size(); i++)
        {
            for(int j = i+1; j < clusters.size(); j++)
            {
                cout << endl;
                cout << "Cluster " << i << ": " << clusters[i].first << " "; 
                cout << "Cluster " << j << ": " << clusters[j].first << "\n"; 
            }
        }
        
        minCoords = searchMin(distanceMatrix);
        distance = (distanceMatrix[minCoords.x][minCoords.y]);
        cout << "Merging Clusters: " << clusters[minCoords.x].first << " and " << clusters[minCoords.y].first << " with distance " << distance << endl;

        copy = distanceMatrix;
        for(int i = 0; i < copy.size(); i++)
        {
                int s1 = clusters[minCoords.x].first.size();
                int s2 = clusters[minCoords.y].first.size();
                int div = s1 + s2;
                distance = (distanceMatrix[i][minCoords.x]*s1 + distanceMatrix[i][minCoords.y]*s2)/div;
                
                if(i!=minCoords.x)
                {
                    copy[i][minCoords.x] = distance;
                    copy[minCoords.x][i] = distance;
                }
        }
        
        copy.erase(copy.begin()+minCoords.y);
        
        for(int i = 0; i < copy.size(); i++)
        {
            copy[i].erase(copy[i].begin()+minCoords.y);
        }
        
        distanceMatrix = copy;
        clusters = reduceClusters(minCoords, clusters);
    }
    cout << "\nNewick Format: " << clusters[0].second;

}

fileContents readFile(string inFile)
{
    fileContents fileValues;
    string line;
    string delim = "\t";
    string token;
    size_t pos = 0;
    bool readSize = false;
    bool readNodes = false;
    int scoresLineCounter = 0;
    int scoresIndexCounter = 0;
    string::size_type sz;

    ifstream file;
    file.open(inFile);
    if(file.is_open())
    {
        while(getline(file, line))
        {
            if(!readSize) //if description has not been read
            {
                fileValues.size = stoi(line); 
                readSize = true;
            }
            else if(!readNodes) //if AASequence has not been read
            {
                while((pos = line.find(delim)) != string::npos)
                {
                    token = line.substr(0,pos);
                    fileValues.nodes.push_back(token);
                    line.erase(0, pos + delim.length());
                }
                fileValues.nodes.push_back(line);
                readNodes = true;
            }
            else //Read Matrix Scores
            {
                scoresIndexCounter = 0;
                while((pos = line.find(delim)) != string::npos)
                {
                    
                    token = line.substr(0,pos);
                    if(scoresLineCounter == 0)                                
                    {
                        vector<double> score;
                        score.push_back(stod(token, &sz));
                        fileValues.matrix.push_back(score);
                    }
                    else
                    {
                        fileValues.matrix[scoresIndexCounter].push_back(stod(token, &sz));
                    }
                    
                    line.erase(0, pos + delim.length());
                    scoresIndexCounter++;
                }
                if(scoresLineCounter == 0)
                {
                        vector<double> score;
                        score.push_back(stod(line, &sz));
                        fileValues.matrix.push_back(score);
                }
                else
                {
                        fileValues.matrix[scoresIndexCounter].push_back(stod(line, &sz));
                }
                scoresLineCounter++;
            }
        }

        file.close();
    }
    else
    {
        cerr << "Error, unable to open file: " << inFile << endl;
    }
    return fileValues;
}

void printVector(vector<auto> data)
{
    for(int i=0; i<data.size(); i++)
    {
        cout << data[i] << " ";
    }
    cout << endl;
}

void printMatrix(vector<vector<auto>> data)
{
    for(int i = 0; i < data.size(); i++)
    {
        for(int j = 0; j < data[i].size(); j++)
        {
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
}

string neighborJoin(fileContents data)
{
    vector<pair<string,string>> clusters = setClusters(data.nodes);
    vector<double> primeValues;
    vector<vector<double>> primeMatrix;
    vector<vector<double>> distanceMatrix = data.matrix;
    coordinates minCoords;
    printMatrix(distanceMatrix);
    while(clusters.size() > 2)
    {

        primeValues = setPrimeValues(distanceMatrix);
        printPrimeValues(clusters, primeValues);
        primeMatrix = setPrimeMatrix(distanceMatrix, primeValues);
        printPrimeMatrix(clusters, primeMatrix);
        minCoords = searchMin(primeMatrix);      
        cout << "Merging Clusters: " << clusters[minCoords.x].first << " and " << clusters[minCoords.y].first << endl;
        double xdist = (0.5*distanceMatrix[minCoords.x][minCoords.y])+(0.5*(primeValues[minCoords.x] - primeValues[minCoords.y]));
        double ydist = (0.5*distanceMatrix[minCoords.x][minCoords.y])+(0.5*(primeValues[minCoords.y] - primeValues[minCoords.x]));
        cout << "\tDistance between " << clusters[minCoords.x].first << " and ancestral node = " << xdist << endl;
        cout << "\tDistance between " << clusters[minCoords.y].first << " and ancestral node = " << ydist << endl;
        clusters = reduceClusters(minCoords, clusters);
        distanceMatrix = reduceDistanceMatrix(minCoords, distanceMatrix, "n");
    }

    minCoords.x = 0;
    minCoords.y = 1;
    cout << "Merging Clusters: " << clusters[minCoords.x].first << " and " << clusters[minCoords.y].first << endl;
    double xdist = (0.5*distanceMatrix[minCoords.x][minCoords.y])+(0.5*(primeValues[minCoords.x] - primeValues[minCoords.y]));
    double ydist = (0.5*distanceMatrix[minCoords.x][minCoords.y])+(0.5*(primeValues[minCoords.y] - primeValues[minCoords.x]));
    cout << "\tDistance between " << clusters[minCoords.x].first << " and ancestral node = " << xdist << endl;
    cout << "\tDistance between " << clusters[minCoords.y].first << " and ancestral node = " << ydist << endl;

    clusters = reduceClusters(minCoords, clusters);

    cout << clusters[0].second; 
}

coordinates searchMin(vector<vector<double>> matrix)
{
    coordinates coords;
    coords.x = 0;
    coords.y = 1;

    double smallest = matrix[0][1];
    double nextValue;

    for(int i = 0; i < matrix.size(); i++)
    {
        for(int j = i+1; j < matrix[0].size(); j++)
        {
            nextValue = matrix[i][j];
            if(nextValue < smallest)
            {
                smallest = nextValue;
                coords.x = i;
                coords.y = j;
            }   
        }
    }
    return coords;
}

vector<double> setPrimeValues(vector<vector<double>> matrix)
{
    vector<double> primes;
    double rowSum;
    double n;

    for(int i = 0; i < matrix.size(); i++)
    {
        rowSum = 0;
        for(int j = 0; j < matrix[0].size(); j++)
        {
            rowSum += matrix[i][j];
        }
        n = matrix.size() - 2;
        primes.push_back(rowSum/n);
    }
    return primes;
}

vector<vector<double>> setPrimeMatrix(vector<vector<double>> matrix, vector<double> primes)
{
    vector<vector<double>> primeMatrix = matrix;
    
    for(int i = 0; i < matrix.size(); i++)
    {
        for(int j = 0; j < matrix[0].size(); j++)
        {
            if(i!=j)
                primeMatrix[i][j] = matrix[i][j] - primes[i] - primes[j];

        }
    }
    return primeMatrix;
}

vector<vector<double>> reduceDistanceMatrix(coordinates coords, vector<vector<double>> matrix, string algo)
{
    vector<vector<double>> copy = matrix;
    double distance;

   for(int i = 0; i < copy.size(); i++)
   {
        if(algo == "n")
            distance = (matrix[i][coords.x] + matrix[i][coords.y] - matrix[coords.x][coords.y])/2;
        else
            distance = (matrix[i][coords.x] + matrix[i][coords.y])/2;

        //cout << "dist: " << distance << endl;
        if(i!=coords.x)
        {
            copy[i][coords.x] = distance;
            copy[coords.x][i] = distance;
        }
   }
   
    copy.erase(copy.begin()+coords.y);
    
    for(int i = 0; i < copy.size(); i++)
    {
        copy[i].erase(copy[i].begin()+coords.y);
    }

    return copy;
}

vector<pair<string,string>> setClusters(vector<string> data)
{
    vector<pair<string,string>> clusters;
    pair<string,string> next;
    for(int i = 0; i < data.size(); i++)
    {
        next.first = data[i];
        next.second = data[i];
        clusters.push_back(next);
    }
    return clusters;
}

vector<pair<string,string>> reduceClusters(coordinates coords, vector<pair<string,string>> clusters)
{
    string newick;
    newick += "(";
    newick += clusters[coords.x].second;
    newick += ",";
    newick += clusters[coords.y].second;
    newick += ")";
    clusters[coords.x].first+=clusters[coords.y].first;
    clusters[coords.x].second = newick;
    clusters.erase(clusters.begin() + coords.y);
    return clusters;
}

void printPrimeValues(vector<pair<string,string>> clusters, vector<double> primes)
{
    cout << "R values" << endl;
    for(int i = 0; i < clusters.size(); i++)
    {
        cout << clusters[i].first << " " << primes[i] << " ";
    }
    cout << endl;
}

void printPrimeMatrix(vector<pair<string,string>> clusters, vector<vector<double>> primeMatrix)
{
    cout << "TD matrix:\n";
    for(int i = 0; i < primeMatrix.size(); i++)
    {
        for(int j = i+1; j < primeMatrix.size(); j++)
        {
            cout << "TD(" << clusters[i].first << "," << clusters[j].first << ")";
            cout << " = " << primeMatrix[i][j] << ",\t";
        }
        cout << endl;
    }
}