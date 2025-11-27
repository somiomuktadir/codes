#ifndef MENU_H
#define MENU_H

namespace Menu {
    void showMainMenu();
    
    // Sub-menus
    void showCombinatoricsMenu();
    void showProbabilityMenu();
    void showDistributionsMenu();
    
    // Helpers
    void clearScreen();
    void pause();
}

#endif // MENU_H
