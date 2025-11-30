#include "TestFramework.h"

int main() {
    return TestFramework::TestRegistry::getInstance().runAll();
}
