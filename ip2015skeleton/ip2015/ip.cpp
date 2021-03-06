#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <stack>


bool comparePixel(Pixel &pix1, Pixel &pix2){
    bool similar = true;
    for (int i = 0; i < 3; ++i) similar &= (abs(pix1.getColor(i) - pix2.getColor(i)) < .05);
    return similar;
}


void checkNeighbors(int x, int y, Image &src, vector<bool>* neighbors){
    int srcWidth = src.getWidth();
    int srcHeight = src.getHeight();
    
    Pixel pix = src.getPixel(x, y);
    Pixel pixN;
    
    if(x > 0){
        if (y > 0) {
            pixN = src.getPixel(x-1,y-1);
            neighbors->at(0) = comparePixel(pix,pixN);
        }
        pixN = src.getPixel(x-1,y);
        neighbors->at(3) = comparePixel(pix, pixN);
        if (y < srcHeight - 1) {
            pixN = src.getPixel(x-1,y+1);
            neighbors->at(5) = comparePixel(pix, pixN);
        }
    }
    if (x < srcWidth - 1) {
        if (y > 0) neighbors->at(2) = comparePixel(pix, src.getPixel(x+1, y-1, pixN));
        neighbors->at(4) = comparePixel(pix, src.getPixel(x+1, y, pixN));
        if (y < srcHeight - 1) neighbors->at(7) = comparePixel(pix, src.getPixel(x+1, y+1, pixN));
    }
    if(y > 0) neighbors->at(1) = comparePixel(pix, src.getPixel(x, y-1, pixN));
    if(y < srcHeight -1) neighbors->at(6) = comparePixel(pix, src.getPixel(x, y+1, pixN));
     
    
}

void drawLine(int startX, int startY, int endX, int endY, int blockSize, Image &src){
    bool xChanges = startX != endX;
    bool yChanges = startY != endY;
    Pixel red(1,0,0);
    if(xChanges){
        if(yChanges){
            if(endX > startX){
                if(endY > startY){
                    for(int i = 0; i < endX - startX; i++){
                        src.setPixel(startX+i, startY+i, red);
                    }
                }else{
                    for(int i = 0; i<endX - startX; i++){
                        src.setPixel(startX+i, startY-i, red);
                    }
                }
            }else{
                if(endY > startY){
                    for(int i = 0; i < endX - startX; i++){
                        src.setPixel(startX-i, startY+i, red);
                    }
                }else{
                    for(int i = 0; i<endX - startX; i++){
                        src.setPixel(startX-i, startY-i, red);
                    }
                }
            }
        }else{
            int i = min(startX, endX);
            int bigI = max(startX, endX);
            int j = startY;
            while(i <= bigI){
                src.setPixel(i, j, red);
                i++;
            }
        }
    }else if(yChanges){
        int j = min(startY, endY);
        int bigJ = max(startY, endY);
        int i = startX;
        while(j <= bigJ){
            src.setPixel(i, j, red);
            j++;
        }
    }
}



void drawEdge(int startX, int startY, int endX, int endY, int blockSize, Image &src){
    int startCenterX = startX*blockSize + blockSize/2;
    int startCenterY = startY*blockSize + blockSize/2;
    int endCenterY = endY*blockSize + blockSize/2;
    int endCenterX = endX*blockSize + blockSize/2;
    drawLine(startCenterX,startCenterY,endCenterX,endCenterY,blockSize,src);
}


void getNeighbor(int startX, int startY, int* neighborX, int* neighborY, int index){
    switch(index){
        case 0: {
            *neighborX = startX-1;
            *neighborY = startY-1;
            break;
        }
        case 1: {
            *neighborX = startX;
            *neighborY = startY-1;
            break;
        }
        case 2: {
            *neighborX = startX+1;
            *neighborY = startY-1;
            break;
        }
        case 3: {
            *neighborX = startX-1;
            *neighborY = startY;
            break;
        }
        case 4: {
            *neighborX = startX+1;
            *neighborY = startY;
            break;
        }
        case 5: {
            *neighborX = startX-1;
            *neighborY = startY+1;
            break;
        }
        case 6: {
            *neighborX = startX;
            *neighborY = startY+1;
            break;
        }
        case 7: {
            *neighborX = startX+1;
            *neighborY = startY+1;
            break;
        }
    }
}

void getClockWiseNeighbor(int startX, int startY, int* neighborX, int* neighborY, int index){
    switch(index){
        case 0: {
            *neighborX = startX-1;
            *neighborY = startY-1;
            break;
        }
        case 1: {
            *neighborX = startX;
            *neighborY = startY-1;
            break;
        }
        case 2: {
            *neighborX = startX+1;
            *neighborY = startY-1;
            break;
        }
        case 3: {
            *neighborX = startX+1;
            *neighborY = startY;
            break;
        }
        case 4: {
            *neighborX = startX+1;
            *neighborY = startY+1;
            break;
        }
        case 5: {
            *neighborX = startX;
            *neighborY = startY+1;
            break;
        }
        case 6: {
            *neighborX = startX-1;
            *neighborY = startY+1;
            break;
        }
        case 7: {
            *neighborX = startX-1;
            *neighborY = startY;
            break;
        }
    }
}


int getNeighborIndex(int startIndex, int blockSizeX, int blockSizeY, int neighborIndex){
    int resultIndex = -1;
    switch(neighborIndex){
        case 0: {
            resultIndex = startIndex - blockSizeY - 1;
            break;
        }
        case 1: {
            resultIndex = startIndex - 1;
            break;
        }
        case 2: {
            resultIndex = startIndex + blockSizeY - 1;
            break;
        }
        case 3: {
            resultIndex = startIndex - blockSizeY;
            break;
        }
        case 4: {
            resultIndex = startIndex + blockSizeY;
            break;
        }
        case 5: {
            resultIndex = startIndex - blockSizeY + 1;
            break;
        }
        case 6: {
            resultIndex = startIndex + 1;
            break;
        }
        case 7: {
            resultIndex = startIndex + blockSizeY + 1;
            break;
        }
    }
    return resultIndex;
}

// Pixel at i,j represents left corner of 2x2 block to check
void checkRedundancy(int i,
                     int j,
                     int blockSizeX,
                     int blockSizeY,
                     vector< vector<bool> > * similarityGraph) {
    vector<bool> adjList = similarityGraph->at(i*blockSizeY + j);
    bool upperLeft = adjList[4] && adjList[6] && adjList[7];
    adjList = similarityGraph->at(i*blockSizeY + j + 1);
    bool bottomLeft = adjList[2] && adjList[4];
    adjList = similarityGraph->at((i+1)*blockSizeY + j + 1);
    bool bottomRight = adjList[1];
    if (upperLeft && bottomRight && bottomLeft) {
        similarityGraph->at(i*blockSizeY + j)[7] = false;
        similarityGraph->at((i+1)*blockSizeY + j+1)[0] = false;
        similarityGraph->at(i*blockSizeY + j+1)[2] = false;
        similarityGraph->at((i+1)*blockSizeY + j)[5] = false;
    }
}

/*
 * Borders is: 
 * {left, right, top, bottom}
 */
bool bounds(int i,
            int j,
            int * borders) {
    return (borders[0] <= i && i <= borders[1] && borders[2] <= j && j <= borders[3]);
}

/*
 * Returns the size of the connected component at pixel i, j.
 */
int sizeConnectedComponent(int i,
                           int j,
                           int blockSizeX,
                           int blockSizeY,
                           int * borders,
                           vector< vector<bool> > * similarityGraph) {
    vector<bool> visited(blockSizeX*blockSizeY, false);
    stack<int> nodesToVisit;
    vector<bool> adjlist = similarityGraph->at(i*blockSizeY+j);
    int neighborX = i;
    int neighborY = j;
    int numberOfConnectedComponents = 0;
    for (int index = 0; index < adjlist.size(); ++index) {
        if(adjlist[index] && bounds(i, j, borders)) {
            getNeighbor(i, j, &neighborX, &neighborY , index);
            nodesToVisit.push(neighborX * blockSizeY + neighborY);
            visited[neighborX * blockSizeY + neighborY] = true;
            numberOfConnectedComponents++;
        }
    }

    int neighbor = i * blockSizeY + blockSizeX;
    while (!nodesToVisit.empty()) {
        int visitedNode = nodesToVisit.top();
        adjlist = similarityGraph->at(visitedNode);
        for (int index = 0; index < adjlist.size(); ++index) {
            getNeighbor(visitedNode / blockSizeY, visitedNode % blockSizeY, &neighborX, &neighborY, index);
            neighbor = neighborX * blockSizeY + neighborY;
            if (adjlist[index] && !visited[neighbor] && bounds(neighborX, neighborY, borders)) {
                nodesToVisit.push(neighbor);
                visited[neighbor] = true;
                numberOfConnectedComponents++;
            }
            
        }
        nodesToVisit.pop();
    }
    
    return numberOfConnectedComponents;
}

void sparsePixels(int i,
                  int j,
                  int blockSizeX,
                  int blockSizeY,
                  int * weightMajor,
                  int * weightMinor,
                  vector< vector<bool> > * similarityGraph) {
    int left = max(0, i - 3);
    int right = min(blockSizeX, i + 4);
    int top = max(0, j - 3);
    int bottom = min(blockSizeY, j + 4);
    int borders[] = {left, right, top, bottom};
    int major = sizeConnectedComponent(i + 1, j, blockSizeX, blockSizeY, borders, similarityGraph);
    int minor = sizeConnectedComponent(i, j, blockSizeX, blockSizeY, borders, similarityGraph);
    
    if (major > minor) *weightMinor += major - minor;
    else *weightMajor += minor - major;
}

int getNodeDegree(vector<bool> & neighbors) {
    int total = 0;
    for (int i = 0; i < neighbors.size(); ++i) {
        if (neighbors[i]) ++total;
    }
    return total;
}

/*
 * Check if one of the diagonal nodes are of valence one. If so, add 5 to
 * its weight.
 */
void islands(int i,
             int j,
             int blockSizeX,
             int blockSizeY,
             int * weightMajor,
             int * weightMinor,
             vector< vector<bool> > & similarityGraph) {
    int topLeft = getNodeDegree(similarityGraph[i*blockSizeY+j]);
    int topRight = getNodeDegree(similarityGraph[i*blockSizeY+j+1]);
    int botLeft = getNodeDegree(similarityGraph[(i+1)*blockSizeY+j]);
    int botRight =getNodeDegree(similarityGraph[(i+1)*blockSizeY+j+1]);
    if (topLeft == 1 || botRight == 1) weightMinor += 5;
    //question
    if (topRight == 1 || botLeft == 1) weightMajor += 5;
}

/*
 *goes to the next node given a node index and a visited path. Presuming that the node in question is a valence-2 node
 */
int goToNextNode(int nodeIndex, int blockSizeX, int blockSizeY, int visited, vector< vector<bool> > & similarityGraph){
    vector<bool> neighbors = similarityGraph[nodeIndex];
    for(int i = 0; i < neighbors.size(); i++){
        if(neighbors[i]){
            int neighborIndex = getNeighborIndex(nodeIndex, blockSizeX, blockSizeY, i);
            if(neighborIndex != visited){
                return neighborIndex;
            }
        }
    }
    return -1;
}

void curve(int i,
           int j,
           int blockSizeX,
           int blockSizeY,
           int * weightMajor,
           int * weightMinor,
           vector< vector<bool> > & similarityGraph) {
    int topLeftIndex = i*blockSizeY+j;
    int topRightIndex = i*blockSizeY+j+1;
    int botLeftIndex = (i+1)*blockSizeY+j;
    int botRightIndex = (i+1)*blockSizeY+j+1;
    int topLeft = getNodeDegree(similarityGraph[topLeftIndex]);
    int topRight = getNodeDegree(similarityGraph[topRightIndex]);
    int botLeft = getNodeDegree(similarityGraph[botLeftIndex]);
    int botRight =getNodeDegree(similarityGraph[botRightIndex]);
    int topLeftVisited = botRightIndex;
    int topRightVisited = botLeftIndex;
    int botLeftVisited = topRightIndex;
    int botRightVisited = topLeftIndex;
    int newTopLeftIndex, newTopRightIndex, newBotLeftIndex, newBotRightIndex;
    int majorLength = 1;
    int minorLength = 1;
    while(topLeft == 2){
        minorLength++;
        newTopLeftIndex = goToNextNode(topLeftIndex, blockSizeX, blockSizeY, topLeftVisited, similarityGraph);
        topLeftVisited = topLeftIndex;
        topLeftIndex = newTopLeftIndex;
        topLeft = getNodeDegree(similarityGraph[newTopLeftIndex]);
    }
    while(botRight == 2){
        minorLength++;
        newBotRightIndex = goToNextNode(botRightIndex, blockSizeX, blockSizeY, botRightVisited, similarityGraph);
        botRightVisited = botRightIndex;
        botRightIndex = newBotRightIndex;
        botRight = getNodeDegree(similarityGraph[newBotRightIndex]);
    }
    while(topRight == 2){
        majorLength++;
        newTopRightIndex = goToNextNode(topRightIndex, blockSizeX, blockSizeY, topRightVisited, similarityGraph);
        topRightVisited = topRightIndex;
        topRightIndex = newTopRightIndex;
        topRight = getNodeDegree(similarityGraph[newTopRightIndex]);
    }
    while(botLeft == 2){
        majorLength++;
        newBotLeftIndex = goToNextNode(botLeftIndex, blockSizeX, blockSizeY, botLeftVisited, similarityGraph);
        botLeftVisited = botLeftIndex;
        botLeftIndex = newBotLeftIndex;
        botLeft = getNodeDegree(similarityGraph[newBotLeftIndex]);
    }
    if(majorLength>minorLength){
        weightMajor += 6;
    }else if(minorLength>majorLength){
        weightMinor += 6;
    }
}

/*
 * Pixel at i,j is the left hand corner of a 2x2 block with diagonals.
 */
void chooseDiagonals(int i,
                     int j,
                     int blockSizeX,
                     int blockSizeY,
                     vector< vector<bool> > * similarityGraph) {
    int weightMajor = 0;
    int weightMinor = 0;
    sparsePixels(i, j, blockSizeX, blockSizeY, &weightMajor, &weightMinor, similarityGraph);
    islands(i, j, blockSizeX, blockSizeY, &weightMajor, &weightMinor, *similarityGraph);
    curve(i, j, blockSizeX, blockSizeY, &weightMajor, &weightMinor, *similarityGraph);
    if (weightMinor > weightMajor) {
        similarityGraph->at((i)*blockSizeY+j+1)[2] = false;
        similarityGraph->at((i+1)*blockSizeY+j)[5] = false;
    } else {
        similarityGraph->at((i)*blockSizeY+j)[7] = false;
        similarityGraph->at((i+1)*blockSizeY+j+1)[0] = false;
    }
}

void drawLineFromMidpoint(int i, int j, int pixelSize, int neighborIndex, Image& src, Pixel p){
    int halfRange = pixelSize/4;
    if(neighborIndex == 0 || neighborIndex == 7){
        for(int k = - halfRange; k <= halfRange; k++){
            src.setPixel(i + k, j - k, p);
            src.setPixel(i + k + 1, j - k, p);
            src.setPixel(i + k - 1, j - k, p);
        }
    }else if(neighborIndex == 2 || neighborIndex == 5){
        for(int k = - halfRange ; k <= halfRange; k++){
            src.setPixel(i + k, j + k, p);
            src.setPixel(i + k - 1, j + k, p);
            src.setPixel(i + k + 1, j + k, p);
        }
    }
}

void reshapePixel(int i, int j, int pixelSize, int neighborIndex, Image& src){
    
    Pixel p = src.getPixel(i*pixelSize + pixelSize/2+1, j*pixelSize + pixelSize/2+2);
    if(neighborIndex == 0){
        int startDiagX = i*pixelSize;
        int startDiagY = j*pixelSize;
        for(int k = 0; k <= pixelSize/4+1; k++){
            drawLineFromMidpoint(startDiagX+k, startDiagY+k, pixelSize, 0, src, p);
        }
    }else if(neighborIndex == 2){
        int startDiagX = i*pixelSize + pixelSize;
        int startDiagY = j*pixelSize;
        for(int k = 0; k <= pixelSize/4+1; k++){
            drawLineFromMidpoint(startDiagX-k, startDiagY+k, pixelSize, 2, src, p);
        }
    }else if(neighborIndex == 5){
        int startDiagX = i*pixelSize;
        int startDiagY = j*pixelSize + pixelSize;
        for(int k = 0; k <= pixelSize/4+1; k++){
            drawLineFromMidpoint(startDiagX+k, startDiagY-k, pixelSize, 5, src, p);
        }
    }else if(neighborIndex == 7){
        int startDiagX = i*pixelSize + pixelSize;
        int startDiagY = j*pixelSize + pixelSize;
        for(int k = 0; k <= pixelSize/4+1; k++){
            drawLineFromMidpoint(startDiagX-k, startDiagY-k, pixelSize, 7, src, p);
        }
    }
}

void reshapePixels(vector< vector<bool> > &similarityGraph,
                 Image& src,
                 int blockSizeX,
                 int blockSizeY,
                 int pixelSize){
    for(int i = 0; i < blockSizeX; i++){
        for(int j = 0; j < blockSizeY; j++){
            vector<bool> adjlist = similarityGraph[i*blockSizeY+j];
            for (int k = 0; k < adjlist.size(); ++k)
                if (adjlist[k]) reshapePixel(i, j, pixelSize, k, src);
        }
    }
}

void reshapePixelGl(int i,
                    int j,
                    vector<bool> adjList,
                    ImageGraph* src,
                    int blockSizeX,
                    int blockSizeY,
                    int pixelSize){

    ImagePixel* p = &currentImageGraph->at(i*blockSizeY+j);

    ImagePixel* T;
    ImagePixel* TL;
    ImagePixel* TR;
    ImagePixel* L;
    ImagePixel* R;
    
    if(adjList[0]){
        TL = &(*src)[(i-1)*blockSizeY + (j-1)];
        T = &(*src)[i*blockSizeY + (j-1)];
        L = &(*src)[(i-1)*blockSizeY + j];

        p->vertices[0][0] += pixelSize/4;
        p->vertices[0][1] += pixelSize/4;
        p->vertices[1][0] -= pixelSize/4;
        p->vertices[1][1] -= pixelSize/4;
        
        
        T->vertices[2][0] += 0.25*pixelSize;
        T->vertices[2][1] += 0.25*pixelSize;
        T->vertices[3][0] += 0.25*pixelSize;
        T->vertices[3][1] += 0.25*pixelSize;
        
        TL->vertices[4][0] -= 0.25*pixelSize;
        TL->vertices[4][1] -= 0.25*pixelSize;
        TL->vertices[5][0] += 0.25*pixelSize;
        TL->vertices[5][1] += 0.25*pixelSize;
        
        L->vertices[6][0] -= 0.25*pixelSize;
        L->vertices[6][1] -= 0.25*pixelSize;
        L->vertices[7][0] -= 0.25*pixelSize;
        L->vertices[7][1] -= 0.25*pixelSize;
    }
    
    if(adjList[2]){
        T = &(*src)[i*blockSizeY + (j-1)];
        TR = &(*src)[(i+1)*blockSizeY + (j-1)];
        R = &(*src)[(i+1)*blockSizeY + j];
        
        p->vertices[6][0] += 0.25*pixelSize;
        p->vertices[6][1] -= 0.25*pixelSize;
        p->vertices[7][0] -= 0.25*pixelSize;
        p->vertices[7][1] += 0.25*pixelSize;

        T->vertices[4][0] -= 0.25*pixelSize;
        T->vertices[4][1] += 0.25*pixelSize;
        T->vertices[5][0] -= 0.25*pixelSize;
        T->vertices[5][1] += 0.25*pixelSize;
        
        TR->vertices[2][0] -= 0.25*pixelSize;
        TR->vertices[2][1] += 0.25*pixelSize;
        TR->vertices[3][0] += 0.25*pixelSize;
        TR->vertices[3][1] -= 0.25*pixelSize;
        
        R->vertices[0][0] += 0.25*pixelSize;
        R->vertices[0][1] -= 0.25*pixelSize;
        R->vertices[1][0] += 0.25*pixelSize;
        R->vertices[1][1] -= 0.25*pixelSize;
    }
    
    
}

void reshapePixelsGl(vector< vector<bool> > &similarityGraph,
                     ImageGraph* src,
                     int blockSizeX,
                     int blockSizeY,
                     int pixelSize){
    for(int i=0; i<blockSizeX; i++){
        for(int j=0; j<blockSizeY; j++){
            reshapePixelGl(i, j, similarityGraph[i*blockSizeY+j], src, blockSizeX, blockSizeY, pixelSize);
        }
    }

}

vector<vector<bool>> simplifiedSimilarityGraph(vector<vector<bool>> &similarityGraph,
                                               int blockSizeX,
                                               int blockSizeY){
    vector<vector<bool>> simplifiedGraph = similarityGraph;
    for(int i = 0; i<blockSizeX; i++){
        for(int j = 0; j<blockSizeY; j++){
            if(simplifiedGraph[i*blockSizeY+j][1]){
                if(similarityGraph[i*blockSizeY+j][3] && similarityGraph[i*blockSizeY+j][4]){
                    simplifiedGraph[i*blockSizeY+j][1] = false;
                    simplifiedGraph[i*blockSizeY+j-1][6] = false;
                }else if(similarityGraph[i*blockSizeY+j-1][3] && similarityGraph[i*blockSizeY+j-1][4]){
                    simplifiedGraph[i*blockSizeY+j][1] = false;
                    simplifiedGraph[i*blockSizeY+j-1][6] = false;
                }
            }
            if(simplifiedGraph[i*blockSizeY+j][3]){
                if(similarityGraph[i*blockSizeY+j][1] && similarityGraph[i*blockSizeY+j][6]){
                    simplifiedGraph[i*blockSizeY+j][3] = false;
                    simplifiedGraph[(i-1)*blockSizeY+j][4] = false;
                }else if(similarityGraph[(i-1)*blockSizeY+j][1] && similarityGraph[(i-1)*blockSizeY+j][6]){
                    simplifiedGraph[i*blockSizeY+j][3] = false;
                    simplifiedGraph[(i-1)*blockSizeY+j][4] = false;
                }
            }
        }
    }
    return simplifiedGraph;
}


int getClockWiseNodeDegree(vector<bool> & neighbors) {
    int total = 0;
    for (int i = 0; i < 8; ++i) {
        if (neighbors[i]) ++total;
    }
    return total;
}


vector< vector<bool>> clockWiseSimplifiedGraph(vector<vector<bool>> &simplifiedGraph, int blockSizeX, int blockSizeY){
    vector<vector<bool>> clockWiseGraph = *new vector<vector<bool>>(blockSizeX * blockSizeY, vector<bool>(9,false));
    for(int i = 0; i < blockSizeX; i++){
        for(int j = 0; j < blockSizeY; j++){
            clockWiseGraph[i*blockSizeY+j][0] = simplifiedGraph[i*blockSizeY+j][0];
            clockWiseGraph[i*blockSizeY+j][1] = simplifiedGraph[i*blockSizeY+j][1];
            clockWiseGraph[i*blockSizeY+j][2] = simplifiedGraph[i*blockSizeY+j][2];
            clockWiseGraph[i*blockSizeY+j][3] = simplifiedGraph[i*blockSizeY+j][4];
            clockWiseGraph[i*blockSizeY+j][4] = simplifiedGraph[i*blockSizeY+j][7];
            clockWiseGraph[i*blockSizeY+j][5] = simplifiedGraph[i*blockSizeY+j][6];
            clockWiseGraph[i*blockSizeY+j][6] = simplifiedGraph[i*blockSizeY+j][5];
            clockWiseGraph[i*blockSizeY+j][7] = simplifiedGraph[i*blockSizeY+j][3];
            clockWiseGraph[i*blockSizeY+j][8] = false;
        }
    }
    return clockWiseGraph;
}

vector<vector<GLfloat>> extractClockWiseHalf(int i,
                                         int j,
                                         int blockSizeX,
                                         int blockSizeY,
                                         int index){
    ImagePixel* p = &currentImageGraph->at(i*blockSizeY+j);
    
    vector<vector<GLfloat>> clockWiseHalf = *new vector<vector<GLfloat>>(4, vector<GLfloat>(2, -1));
    
    switch (index){
        case 0:{
            clockWiseHalf[0] = p->vertices[0];
            clockWiseHalf[1] = p->vertices[7];
            clockWiseHalf[2] = p->vertices[6];
            clockWiseHalf[3] = p->vertices[5];
            break;
        }
        case 1:{
            clockWiseHalf[0] = p->vertices[7];
            clockWiseHalf[1] = p->vertices[6];
            clockWiseHalf[2] = p->vertices[5];
            clockWiseHalf[3] = p->vertices[4];
            break;
        }
        case 2:{
            clockWiseHalf[0] = p->vertices[6];
            clockWiseHalf[1] = p->vertices[5];
            clockWiseHalf[2] = p->vertices[4];
            clockWiseHalf[3] = p->vertices[3];
            break;
        }
        case 3:{
            clockWiseHalf[0] = p->vertices[5];
            clockWiseHalf[1] = p->vertices[4];
            clockWiseHalf[2] = p->vertices[3];
            clockWiseHalf[3] = p->vertices[2];
            break;
        }
        case 4:{
            clockWiseHalf[0] = p->vertices[4];
            clockWiseHalf[1] = p->vertices[3];
            clockWiseHalf[2] = p->vertices[2];
            clockWiseHalf[3] = p->vertices[1];
            break;
        }
        case 5:{
            clockWiseHalf[0] = p->vertices[3];
            clockWiseHalf[1] = p->vertices[2];
            clockWiseHalf[2] = p->vertices[1];
            clockWiseHalf[3] = p->vertices[0];
            break;
        }
        case 6:{
            clockWiseHalf[0] = p->vertices[2];
            clockWiseHalf[1] = p->vertices[1];
            clockWiseHalf[2] = p->vertices[0];
            clockWiseHalf[3] = p->vertices[7];
            break;
        }
        case 7:{
            clockWiseHalf[0] = p->vertices[1];
            clockWiseHalf[1] = p->vertices[0];
            clockWiseHalf[2] = p->vertices[7];
            clockWiseHalf[3] = p->vertices[6];
            break;
        }
        
    }
    
    
    return clockWiseHalf;
}


vector<vector<GLfloat>> extractCounterClockWiseHalf(int i,
                                             int j,
                                             int blockSizeX,
                                             int blockSizeY,
                                             int index){
    ImagePixel* p = &currentImageGraph->at(i*blockSizeY+j);
    
    vector<vector<GLfloat>> counterClockWiseHalf = *new vector<vector<GLfloat>>(4, vector<GLfloat>(2, -1));
    
    switch (index){
        case 0:{
            counterClockWiseHalf[0] = p->vertices[1];
            counterClockWiseHalf[1] = p->vertices[2];
            counterClockWiseHalf[2] = p->vertices[3];
            counterClockWiseHalf[3] = p->vertices[4];
            break;
        }
        case 1:{
            counterClockWiseHalf[0] = p->vertices[0];
            counterClockWiseHalf[1] = p->vertices[1];
            counterClockWiseHalf[2] = p->vertices[2];
            counterClockWiseHalf[3] = p->vertices[3];
            break;
        }
        case 2:{
            counterClockWiseHalf[0] = p->vertices[7];
            counterClockWiseHalf[1] = p->vertices[0];
            counterClockWiseHalf[2] = p->vertices[1];
            counterClockWiseHalf[3] = p->vertices[2];
            break;
        }
        case 3:{
            counterClockWiseHalf[0] = p->vertices[6];
            counterClockWiseHalf[1] = p->vertices[7];
            counterClockWiseHalf[2] = p->vertices[0];
            counterClockWiseHalf[3] = p->vertices[1];
            break;
        }
        case 4:{
            counterClockWiseHalf[0] = p->vertices[5];
            counterClockWiseHalf[1] = p->vertices[6];
            counterClockWiseHalf[2] = p->vertices[7];
            counterClockWiseHalf[3] = p->vertices[0];
            break;
        }
        case 5:{
            counterClockWiseHalf[0] = p->vertices[4];
            counterClockWiseHalf[1] = p->vertices[5];
            counterClockWiseHalf[2] = p->vertices[6];
            counterClockWiseHalf[3] = p->vertices[7];
            break;
        }
        case 6:{
            counterClockWiseHalf[0] = p->vertices[3];
            counterClockWiseHalf[1] = p->vertices[4];
            counterClockWiseHalf[2] = p->vertices[5];
            counterClockWiseHalf[3] = p->vertices[6];
            break;
        }
        case 7:{
            counterClockWiseHalf[0] = p->vertices[2];
            counterClockWiseHalf[1] = p->vertices[3];
            counterClockWiseHalf[2] = p->vertices[4];
            counterClockWiseHalf[3] = p->vertices[5];
            break;
        }
            
    }
    
    
    return counterClockWiseHalf;
}


int getClockWiseNeighborIndex(int index, int i, int j, int blockSizeY, int blockSizeX){
    int newIndex = -1;
    switch(index){
        case 0:{
            newIndex = (i-1)*blockSizeY+(j-1);
            break;
        }
            
        case 1:{
            newIndex = i*blockSizeY + (j-1);
            break;
        }
            
        case 2:{
            newIndex = (i+1)*blockSizeY + (j-1);
            break;
        }
            
        case 3:{
            newIndex = (i+1)*blockSizeY + j;
            break;
        }
            
        case 4:{
            newIndex = (i+1)*blockSizeY + (j+1);
            break;
        }
            
        case 5:{
            newIndex = i*blockSizeY + (j+1);
            break;
        }
            
        case 6:{
            newIndex = (i-1)*blockSizeY + (j+1);
            break;
        }
            
        case 7:{
            newIndex = (i-1)*blockSizeY + j;
            break;
        }
    }
    return newIndex;
}

int getClockWiseNeighborX(int i, int index){
    if(index == 1 || index == 5){
        return i;
    }else if(index == 0 || index == 7 || index == 6){
        return i-1;
    }else{
        return i+1;
    }
}

int getClockWiseNeighborY(int j, int index){
    if(index == 3 || index == 7){
        return j;
    }else if(index == 0 || index == 1 || index == 2){
        return j-1;
    }else{
        return j+1;
    }
}

bool contains(vector<GLfloat> vertex, vector<vector<GLfloat>> array){
    for(int i = 0; i < array.size(); i++){
        if(vertex[0] == array[i][0] && vertex[1] == array[i][1]){
            return true;
        }
    }
    return false;
}

bool startClockWiseVisible(int i, int j, int blockSizeX, int blockSizeY, int index,vector<vector<bool>>* clockWiseGraph){
    
    return true;
//    ImagePixel* p = &currentImageGraph->at(i*blockSizeY+j);
//    
//    int neighborX = getClockWiseNeighborX(i, index+2 % 8);
//    int neighborY = getClockWiseNeighborY(j, index+2 % 8);
//    if(neighborX < 0 || neighborX >= blockSizeX || neighborY < 0 || neighborY >= blockSizeY){
//        return false;
//    }
//    else{
//        int neighborIndex = getClockWiseNeighborIndex((index+2 % 8), i, j, blockSizeY, blockSizeX);
//        ImagePixel* neighbor = &currentImageGraph->at(neighborIndex);
//        Pixel currentPixel = p->getPixel();
//        Pixel neighborPixel = neighbor->getPixel();
//        return !comparePixel(currentPixel, neighborPixel);
//    }
    
}

bool startCounterClockWiseVisible(int i, int j, int blockSizeX, int blockSizeY, int index, vector<vector<bool>>* clockWiseGraph){
    return true;
//    ImagePixel* p = &currentImageGraph->at(i*blockSizeY+j);
//    int neighborX = getClockWiseNeighborX(i, index-2 % 8);
//    int neighborY = getClockWiseNeighborY(j, index-2 % 8);
//    if(neighborX < 0 || neighborX >= blockSizeX || neighborY < 0 || neighborY >= blockSizeY){
//        return false;
//    }else{
//        ImagePixel* neighbor = &currentImageGraph->at(getClockWiseNeighborIndex((index-2 % 8), i, j, blockSizeY, blockSizeX));
//        Pixel currentPixel = p->getPixel();
//        Pixel neighborPixel = neighbor->getPixel();
//        return !comparePixel(currentPixel, neighborPixel);
//
//    }
}



//extract control points recursively
//extremely long code(may be able to use while loop)
void extractCurveControlPoint(int startIndex,
                              int i,
                              int j,
                              int blockSizeX,
                              int blockSizeY,
                              vector<vector<GLfloat>>* curve1,
                              vector<vector<GLfloat>>* curve2,
                              vector<vector<bool>>* clockWiseGraph){
    
    if(i==17&&j==0){
        cout<<"at 17"<<endl;
    }
    
    //get the end index
    int endIndex = -1;
    for(int x = 0; x < 8; x++){
        if(clockWiseGraph->at(i*blockSizeY+j)[x]){
            if(x != startIndex){
                endIndex = x;
                break;
            }
        }
    }
    
    //extract the four halves of vertices
    vector<vector<GLfloat>> startClockWiseHalf = extractClockWiseHalf(i, j, blockSizeX, blockSizeY, startIndex);
    vector<vector<GLfloat>> startCounterClockWiseHalf = extractCounterClockWiseHalf(i, j, blockSizeX, blockSizeY, startIndex);
    vector<vector<GLfloat>> endClockWiseHalf = extractClockWiseHalf(i, j, blockSizeX, blockSizeY, endIndex);
    vector<vector<GLfloat>> endCounterClockWiseHalf = extractCounterClockWiseHalf(i, j, blockSizeX, blockSizeY, endIndex);
    
    
    if(startIndex > endIndex){
        if(startIndex - endIndex > 4){
            if(startClockWiseVisible(i, j, blockSizeX, blockSizeY, startIndex, clockWiseGraph)){
                for(int i = 0; i < 4; i++){
                    if(contains(startClockWiseHalf[i], endCounterClockWiseHalf)){
                        if(!contains(startClockWiseHalf[i], *curve1)){
                            curve1->push_back(startClockWiseHalf[i]);
                        }
                    }
                }
            }
            if(startCounterClockWiseVisible(i, j, blockSizeX, blockSizeY, startIndex, clockWiseGraph)){
                for(int i = 0; i < 4; i++){
                    if(!contains(startCounterClockWiseHalf[i], *curve2)){
                        curve2->push_back(startCounterClockWiseHalf[i]);
                    }
                }
                for(int i = 3; i >= 0; i--){
                    if(!contains(endClockWiseHalf[i], startCounterClockWiseHalf)){
                        if(!contains(endClockWiseHalf[i], *curve2)){
                            curve2->push_back(endClockWiseHalf[i]);
                        }
                    }
                }
            }
            
        }else{
            //put the intersection of the clockWiseHalf of start and the counterClockWiseHalf of end into curve1
            if(startClockWiseVisible(i, j, blockSizeX, blockSizeY, startIndex, clockWiseGraph)){
                for(int i = 0; i < 4; i++){
                    if(!contains(startClockWiseHalf[i], *curve1)){
                        curve1->push_back(startClockWiseHalf[i]);
                    }
                }
                for(int i = 3; i >= 0; i--){
                    if(!contains(endCounterClockWiseHalf[i], startClockWiseHalf)){
                        if(!contains(endCounterClockWiseHalf[i], *curve1)){
                            curve1->push_back(endCounterClockWiseHalf[i]);
                        }
                    }
                }
            }
            
            //put the union of the counterClockWiseHalf of start and the clockWiseHalf of end into curve2
            if(startCounterClockWiseVisible(i, j, blockSizeX, blockSizeY, startIndex, clockWiseGraph)){
                for(int i = 0; i < 4; i++){
                    if(contains(startCounterClockWiseHalf[i], endClockWiseHalf)){
                        if(!contains(startCounterClockWiseHalf[i], *curve2)){
                            curve2->push_back(startCounterClockWiseHalf[i]);
                        }
                    }
                }
            }
            
            
        }
    }else{
        if(endIndex - startIndex > 4){
            if(startCounterClockWiseVisible(i, j, blockSizeX, blockSizeY, startIndex, clockWiseGraph)){
                for(int i = 0; i < 4; i++){
                    if(contains(startCounterClockWiseHalf[i], endClockWiseHalf)){
                        if(!contains(startCounterClockWiseHalf[i], *curve2)){
                            curve2->push_back(startCounterClockWiseHalf[i]);
                        }
                    }
                }
            }
            
            if(startClockWiseVisible(i, j, blockSizeX, blockSizeY, startIndex, clockWiseGraph)){
                for(int i = 0; i < 4; i++){
                    if(!contains(startClockWiseHalf[i], *curve1)){
                        curve1->push_back(startClockWiseHalf[i]);
                    }
                }
                for(int i = 3; i >= 0; i--){
                    if(!contains(endCounterClockWiseHalf[i], startClockWiseHalf)){
                        if(!contains(endCounterClockWiseHalf[i], *curve1)){
                            curve1->push_back(endCounterClockWiseHalf[i]);
                        }
                    }
                }
            }
            
        }else{
            if(startClockWiseVisible(i, j, blockSizeX, blockSizeY, startIndex, clockWiseGraph)){
                for(int i = 0; i < 4; i++){
                    if(contains(startClockWiseHalf[i], endCounterClockWiseHalf)){
                        if(!contains(startClockWiseHalf[i], *curve1)){
                            curve1->push_back(startClockWiseHalf[i]);
                        }
                    }
                }
            }
            
            if(startCounterClockWiseVisible(i, j, blockSizeX, blockSizeY, startIndex, clockWiseGraph)){
                for(int i = 0; i < 4; i++){
                    if(!contains(startCounterClockWiseHalf[i], *curve2)){
                        curve2->push_back(startCounterClockWiseHalf[i]);
                    }
                }
                for(int i = 3; i >= 0; i--){
                    if(!contains(endClockWiseHalf[i], startCounterClockWiseHalf)){
                        if(!contains(endClockWiseHalf[i], *curve2)){
                            curve2->push_back(endClockWiseHalf[i]);
                        }
                    }
                }
            }
            
        }
    }
    //set visited to be true
    clockWiseGraph->at(i*blockSizeY+j)[8] = true;
    //get the next node
    int nextIndex = getClockWiseNeighborIndex(endIndex, i, j, blockSizeY, blockSizeX);
    //if the next node is also of valence 2, recursively extract the control points out of that as well
    if(getClockWiseNodeDegree(clockWiseGraph->at(nextIndex)) == 2){
        int nextStart = -1;
        if(endIndex>3){
            nextStart = endIndex - 4;
        }else{
            nextStart = endIndex + 4;
        }
        int nextI = getClockWiseNeighborX(i, endIndex);
        int nextJ = getClockWiseNeighborY(j, endIndex);
        extractCurveControlPoint(nextStart, nextI, nextJ, blockSizeX, blockSizeY, curve1, curve2, clockWiseGraph);
    }
}


vector<Curve> extractCurveControlPointStartingAt(int i,
                              int j,
                              int blockSizeX,
                              int blockSizeY,
                              vector<vector<bool>>* clockWiseGraph,
                              bool startsWithFirst){
    
    vector<vector<vector<GLfloat>>> curveSet = *new vector<vector<vector<GLfloat>>>();
    vector<Curve>* curves = new vector<Curve>();
    
    vector<vector<GLfloat>> curve1 = *new vector<vector<GLfloat>>();
    vector<vector<GLfloat>> curve2 = *new vector<vector<GLfloat>>();
    Curve* c1 = new Curve();
    Curve* c2 = new Curve();
    int firstIndex = -1;
    int secondIndex = -1;
    bool foundFirst = false;
    //get first and second connection
    for(int x = 0; x < 8; x++){
        if(clockWiseGraph->at(i*blockSizeY+j)[x]){
            if(foundFirst){
                secondIndex = x;
            }else{
                firstIndex = x;
                foundFirst = true;
            }
        }
    }
    
    if(startsWithFirst){
        extractCurveControlPoint(firstIndex, i, j, blockSizeX, blockSizeY, &curve1, &curve2, clockWiseGraph);
    }else{
        extractCurveControlPoint(secondIndex, i, j, blockSizeX, blockSizeY, &curve1, &curve2, clockWiseGraph);
    }
    
    c1->verticies = curve1;
    c2->verticies = curve2;
    
    c1->setColor(currentImageGraph->at(i*blockSizeY+j).color);
    c2->setColor(currentImageGraph->at(i*blockSizeY+j).color);
    
    curveSet.push_back(curve1);
    curveSet.push_back(curve2);
    
    curves->push_back(*c1);
    curves->push_back(*c2);
    
    return *curves;
}

int startOfCurve(int i, int j, int blockSizeX, int blockSizeY, vector<vector<bool>>* clockWiseGraph){
    vector<bool> neighbors = clockWiseGraph->at(i*blockSizeY+j);
    if(getClockWiseNodeDegree(neighbors) == 2){
        bool foundFirst = false;
        int firstIndex = -1;
        int secondIndex = -1;
        for(int x = 0; x<8; x++){
            if(neighbors[x]){
                if(foundFirst){
                    secondIndex = x;
                }else{
                    firstIndex = x;
                    foundFirst = true;
                }
            }
        }
        int firstConnectedIndex = getClockWiseNeighborIndex(firstIndex, i, j, blockSizeY, blockSizeX);
        int secondConnectedIndex = getClockWiseNeighborIndex(secondIndex, i, j, blockSizeY, blockSizeX);
        if(getClockWiseNodeDegree(clockWiseGraph->at(firstConnectedIndex))!= 2){
            return 1;
        }else if(getClockWiseNodeDegree(clockWiseGraph->at(secondConnectedIndex))!= 2){
            return 2;
        }
    }
    return 0;
}

void extractControlPoints(vector<vector<bool>>* simplifiedGraph,
                                                     int blockSizeX,
                                                     int blockSizeY){
    vector<vector<bool>> clockWiseGraph = clockWiseSimplifiedGraph(*simplifiedGraph, blockSizeX, blockSizeY);

    curveVector = new vector<Curve>();
    
    vector<Curve> currentCurves;
    
    for(int i = 0; i<blockSizeX; i++){
        for(int j = 0; j<blockSizeY; j++){
            if(!clockWiseGraph[i*blockSizeY+j][8]){
                if(startOfCurve(i, j, blockSizeX, blockSizeY, &clockWiseGraph)==1){
                    currentCurves = extractCurveControlPointStartingAt(i, j, blockSizeX, blockSizeY, &clockWiseGraph, true);
                    //curves.push_back(currentCurves[0]);
                    //curves.push_back(currentCur
                    curveVector->push_back(currentCurves[0]);
                    curveVector->push_back(currentCurves[1]);
                } else if(startOfCurve(i, j, blockSizeX, blockSizeY, &clockWiseGraph)==2){
                    currentCurves = extractCurveControlPointStartingAt(i, j, blockSizeX, blockSizeY, &clockWiseGraph, false);
                    //curves.push_back(currentCurves[0]);
                    //curves.push_back(currentCurves[1]);
                    curveVector->push_back(currentCurves[0]);
                    curveVector->push_back(currentCurves[1]);
                }
            }
        }
    }
}



/*
 * define your own filter
 * you need to request any input parameters here, not in control.cpp
 */

Image* ip_misc(Image* src,
               ImageGraph* imageGraph,
               vector<vector<bool>>* similarityGraph,
               const int blockSizeX,
               const int blockSizeY,
               const int pixelSize)
{
    
    int srcWidth = src->getWidth();
    int srcHeight = src->getHeight();
    Image* rawGraph = new Image(blockSizeX, blockSizeY, 8);
    currentImageGraph = new ImageGraph(blockSizeX * blockSizeY);
    
    Pixel srcPixel;
    
    int blockWidth = srcWidth / blockSizeX;
    int blockHeight = srcHeight / blockSizeY;
    
    for(int i = 0; i<blockSizeX; i++){
        for(int j=0; j<blockSizeY; j++){
            Pixel outputPixel(0,0,0);
            int blockEndi = blockWidth * i + blockWidth;
            int blockEndj = blockHeight * j + blockHeight;
            blockEndi = min(blockEndi, srcWidth);
            blockEndj = min(blockEndj, srcHeight);
            outputPixel =
              src->getPixel(i * blockWidth + blockWidth / 2, j * blockHeight + blockHeight / 2);
            
            vector<vector<GLfloat>> vertices(8, vector<GLfloat>(2));
            
            // Polygon ordering counter clockwise
            
            for (int k = 0; k < 4; ++k) {
                currentImageGraph->at(i*blockSizeY+j).vertices[2*k][0] = static_cast<float>((i + (k > 1 ? 1: 0)) * pixelSize);
                currentImageGraph->at(i*blockSizeY+j).vertices[2*k][1] = static_cast<float>(16 * pixelSize - (j + (k==1 || k ==2 ? 1:0)) * pixelSize);
                currentImageGraph->at(i*blockSizeY+j).vertices[2*k+1][0] = static_cast<float>((i + (k > 1 ? 1: 0)) * pixelSize);
                currentImageGraph->at(i*blockSizeY+j).vertices[2*k+1][1] = static_cast<float>(16 * pixelSize - (j + (k==1 || k ==2 ? 1:0)) * pixelSize);
            }
            
            currentImageGraph->at(i*blockSizeY+j).setPixel(outputPixel);
            rawGraph->setPixel(i, j, outputPixel);
            
        }
    }
    
    similarityGraph = new vector<vector<bool>>(blockSizeX * blockSizeY, vector<bool>(8,false));
    
    for (int i = 0; i < blockSizeX; ++i) {
        for (int j = 0; j < blockSizeY; ++j) {
            checkNeighbors(i, j, *rawGraph, &similarityGraph->at(i*blockSizeY + j));
        }
    }
    
    
    Image* rawGraphTest = new Image(blockSizeX*pixelSize, blockSizeY*pixelSize, 3);
    
    Pixel rawPixel(0,0,0);
    
    for(int i=0; i<blockSizeX; ++i){
        for(int j=0; j<blockSizeY; ++j){
            rawPixel = rawGraph->getPixel(i, j);
            for(int x = i*pixelSize; x < i*pixelSize+pixelSize; ++x){
                for(int y = j * pixelSize; y < j*pixelSize + pixelSize; ++y){
                    rawGraphTest->setPixel(x, y, rawPixel);
                }
            }
        }
    }

    for (int i = 0; i < blockSizeX - 1; ++i) {
        for (int j = 0; j < blockSizeY - 1; ++j) {
            // If there are cross edges
            if (similarityGraph->at(i*blockSizeY + j)[7] && similarityGraph->at(i*blockSizeY + j + 1)[2]) {
                checkRedundancy(i,j,blockSizeX, blockSizeY, similarityGraph);
            }
            
        }
    }
    for (int i = 0; i < blockSizeX - 1; ++i) {
        for (int j = 0; j < blockSizeY - 1; ++j) {
            if (similarityGraph->at(i*blockSizeY + j)[7] && similarityGraph->at(i*blockSizeY + j + 1)[2]) {
                chooseDiagonals(i,j,blockSizeX, blockSizeY, similarityGraph);
            }
        }
    }
    
    reshapePixelsGl(*similarityGraph, currentImageGraph, blockSizeX, blockSizeY, pixelSize);

    vector<vector<bool>> simplifiedGraph = simplifiedSimilarityGraph(*similarityGraph, blockSizeX, blockSizeY);
    
    extractControlPoints(&simplifiedGraph, blockSizeX, blockSizeY);
    for(int i = 0; i<blockSizeX; ++i){
        for(int j = 0; j<blockSizeY; ++j){
            for(int k = 0; k<8; ++k){
                if(simplifiedGraph[i*blockSizeY+j][k]){
                    int endX = i;
                    int endY = j;

                    getNeighbor(i, j, &endX, &endY, k);
                    drawEdge(i, j, endX, endY, pixelSize, *rawGraphTest);
                }
            }
        }
        
    }
    return rawGraphTest;
}



