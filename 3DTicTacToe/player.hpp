#ifndef _TICTACTOE3D_PLAYER_HPP_
#define _TICTACTOE3D_PLAYER_HPP_

#include "constants.hpp"
#include "deadline.hpp"
#include "move.hpp"
#include "gamestate.hpp"
#include <vector>
#include <map>

namespace TICTACTOE3D
{

class Player
{

public:
    ///perform a move
    ///\param pState the current state of the board
    ///\param pDue time before which we must have returned
    ///\return the next state the board is in after our move

    std::map<std::string, int> m;

    GameState play(const GameState &pState, const Deadline &pDue);
    int alphabeta(const GameState &pState, int const depth, int alpha, int beta, uint8_t player);
};
}

#endif
