#include "system.h"

#include <iostream>

int main() {
	constexpr int L = 50;
	GridType grid = get_randomized_system(L);
	write_png(grid, "test.png");

	return 0;
}
