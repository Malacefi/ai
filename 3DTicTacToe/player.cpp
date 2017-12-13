// Emelie Eriksson & Pontus Brink
#include "player.hpp"
#include <cstdlib>
#include <cmath> 
#include <algorithm>

namespace TICTACTOE3D
{

int keyFunction(int enemyMarkers, int playerMarkers){
	return 10 * enemyMarkers + 1 * playerMarkers;
}

std::map<int, int> points;

void initializePoints(){
	// If the opponent has 1 marker and you 0
	points[keyFunction(1, 0)] = -10;

	// If the player has 1 marker and the opponen 0
	points[keyFunction(0, 1)] = 10;

	// If the opponent has 2 markers and you 0
	points[keyFunction(2, 0)] = -100;

	// If you have 2 markers and the opponent 0
	points[keyFunction(0, 2)] = 80;

	// If the opponent has 2 markers and the player 1 marker
	points[keyFunction(2, 1)] = 100;

	// If the player has 2 markers and the opponent 1 marker
	points[keyFunction(1, 2)] = -90;

	// If the player has 2 markers and the opponent 2 marker
	points[keyFunction(2, 2)] = 0;

	// If the opponent has 3 markers and the player 0
	points[keyFunction(3, 0)] = -1000;

	//If the opponent has 3 markers and the player 1
	points[keyFunction(3, 1)] = 1100;

	// If the player has 3 markers and the opponent 0
	points[keyFunction(0, 3)] = 1000;

	// If the player has 3 markers and the opponent 1
	points[keyFunction(1, 3)] = -990;

	// If the opponent wins
	points[keyFunction(4, 0)] = -10000; 

	// If you win
	points[keyFunction(0, 4)] = 10000; 
}

int eval(const GameState &pState){
	int score = 0;
	
	if(pState.isXWin())
		return 79691999;
	else if (pState.isOWin())
		return -79691999;
	// ONE LOOP 
	for (int l = 0; l<4; ++l){
		for (int j = 0; j<4; ++j){
			int playerRowCounter = 0;
			int oppRowCounter = 0;

			int playerColCounter = 0;
			int oppColCounter = 0;

			int playerLayerCounter = 0;
			int oppLayerCounter = 0;


			for (int i = 0; i<4; ++i){
				if (pState.at(l,j,i)&CELL_X)
					++playerRowCounter;
				else if(pState.at(l,j,i)&CELL_O)
					++oppRowCounter;

				if (pState.at(i,l,j)&CELL_X)
					++playerColCounter;
				else if(pState.at(i,l,j)&CELL_O)
					++oppColCounter;

				if (pState.at(j,i,l)&CELL_X)
					++playerLayerCounter;
				else if(pState.at(j,i,l)&CELL_O)
					++oppLayerCounter;
			}
			score += points[keyFunction(oppRowCounter, playerRowCounter)];

			score += points[keyFunction(oppColCounter, playerColCounter)];

			score += points[keyFunction(oppLayerCounter, playerLayerCounter)];
		}

		int player1Counter = 0;
		int opp1Counter = 0;

		int player2Counter = 0;
		int opp2Counter = 0;

		int player3Counter = 0;
		int opp3Counter = 0;

		int player4Counter = 0;
		int opp4Counter = 0;

		int player5Counter = 0;
		int opp5Counter = 0;

		int player6Counter = 0;
		int opp6Counter = 0;

		int player1DCounter = 0;
		int opp1DCounter = 0;

		int player2DCounter = 0;
		int opp2DCounter = 0;

		int player3DCounter = 0;
		int opp3DCounter = 0;

		int player4DCounter = 0;
		int opp4DCounter = 0;


		// Diagonals
		for (int i = 0; i<4; ++i){
			if (pState.at(l,i,i)&CELL_X)
				++player1Counter;
			else if(pState.at(l,i,i)&CELL_O)
				++opp1Counter;

			if (pState.at(l,i,3-i)&CELL_X)
				++player2Counter;
			else if(pState.at(l,i,3-i)&CELL_O)
				++opp2Counter;

			if (pState.at(i,l,3-i)&CELL_X)
				++player3Counter;
			else if(pState.at(i,l,3-i)&CELL_O)
				++opp3Counter;

			if (pState.at(i,l,i)&CELL_X)
				++player4Counter;
			else if(pState.at(i,l,i)&CELL_O)
				++opp4Counter;

			if (pState.at(i,3-i,l)&CELL_X)
				++player5Counter;
			else if(pState.at(i,3-i,l)&CELL_O)
				++opp5Counter;

			if (pState.at(i,i,l)&CELL_X)
				++player6Counter;
			else if(pState.at(i,i,l)&CELL_O)
				++opp6Counter;

			// Diagonal 1
			if (pState.at(i,i,i)&CELL_X)
				++player1DCounter;
			else if(pState.at(i,i,i)&CELL_O)
				++opp1DCounter;

			// Diagonal 2
			if (pState.at(3-i,i,i)&CELL_X)
				++player2DCounter;
			else if(pState.at(3-i,i,i)&CELL_O)
				++opp2DCounter;
			
			// Diagonal 3
			if (pState.at(i,3-i,i)&CELL_X)
				++player3DCounter;
			else if(pState.at(i,3-i,i)&CELL_O)
				++opp3DCounter;

			// Diagonal 4
			if (pState.at(3-i,3-i,i)&CELL_X)
				++player4DCounter;
			else if(pState.at(3-i,3-i,i)&CELL_O)
				++opp4DCounter;
		}
		score += points[keyFunction(opp1Counter, player1Counter)];
		score += points[keyFunction(opp2Counter, player2Counter)];
		score += points[keyFunction(opp3Counter, player3Counter)];
		score += points[keyFunction(opp4Counter, player4Counter)];
		score += points[keyFunction(opp5Counter, player5Counter)];
		score += points[keyFunction(opp6Counter, player6Counter)];

		score += points[keyFunction(opp1DCounter, player1DCounter)];
		score += points[keyFunction(opp2DCounter, player2DCounter)];
		score += points[keyFunction(opp3DCounter, player3DCounter)];
		score += points[keyFunction(opp4DCounter, player4DCounter)];
	}
	//std::cerr << score << std::endl;

	return score;
}

int sortHelper(GameState &larger, GameState &smaller){
	return eval(larger)>eval(smaller);
}

class Sorter {
	std::map<std::string, int> &m;
public:
	Sorter(std::map<std::string, int> &m) : m(m){}

	bool operator()(GameState &larger, GameState &smaller){
		int big = 0;
		int small = 0;
		if (m.find(larger.toMessage().substr(0,64)) != m.end()){
			big = m[larger.toMessage().substr(0,64)];
		}
		else
			big = eval(larger);
		if (m.find(smaller.toMessage().substr(0,64)) != m.end()){
			small = m[smaller.toMessage().substr(0,64)];
		}
		else
			small = eval(smaller);
		return big > small;
	}
};

int Player::alphabeta(const GameState &pState, int const depth, int alpha, int beta, uint8_t player){

	std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);
    int v = 0;
    //std::cerr << "IN ALPHABETA" << std::endl;
    if (lNextStates.size() == 0 || depth == 0){
    	//std::cerr << "No next states :O"<< std::endl;
    	if (m.find(pState.toMessage().substr(0,64)) != m.end()){
    		//std::cerr << "VALUE IS IN THINGY YAHOOOO" << std::endl;
        	v = m[pState.toMessage().substr(0,64)];
        }
        else{
        	//std::cerr << "VALUE IS NOT IN THINGY TimE TO EVAL" << std::endl;
        	v = eval(pState);
        }
        	
    }
    else if(player == CELL_X){
        v = -10000000;
        for (unsigned i = 0; i<lNextStates.size(); ++i){
            v = std::max(v, alphabeta(lNextStates[i], depth-1,alpha,beta,CELL_O));
            alpha = std::max(alpha, v);
            if (beta <= alpha)
                break;  //Beta prune
        }
    }
    else { //player == CELL_0
        v = 100000000;
        for (unsigned i = 0; i<lNextStates.size(); ++i){
            v = std::min(v, alphabeta(lNextStates[i], depth-1,alpha,beta,CELL_X));
            beta = std::min(beta, v);
            if (beta <= alpha)
                break;  //Alpha prune
        }
    }
    m[pState.toMessage().substr(0,64)] = v;
    return v;
}


GameState Player::play(const GameState &pState,const Deadline &pDue)
{
	initializePoints();
    ////std::cerr << "Processing " << pState.toMessage() << std::endl;

    std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);
    //std::cerr << lNextStates[0].getNextPlayer()<< std::endl;
    std::sort(lNextStates.begin(), lNextStates.end(),Sorter(m));
    //order lNextStates

    if (lNextStates.size() == 0) return GameState(pState, Move());
    int best = 0;
    int bestState = 0;

    for (unsigned i = 0; i<lNextStates.size()/15; ++i){
        int v = alphabeta(lNextStates[i], 0, -10000000, 100000000, CELL_O);
        
    	if (best < v){
            best = v;
            bestState = i;
        }
    }    
    return lNextStates[bestState];
}
}
