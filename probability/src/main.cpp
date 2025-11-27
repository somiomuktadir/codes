#include "Menu.h"
#include <iostream>

int main() {
    try {
        Menu::showMainMenu();
    } catch (const std::exception& e) {
        std::cerr << "Fatal Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
