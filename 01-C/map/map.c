#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Definieren Sie ein enum cardd
typedef enum { N=1, W=2, E=4, S=8 } cardd;

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält
static cardd map[3][3];

char* get_dir(cardd dir) {
    switch (dir) {
        case N: return "N";
        case E: return "E";
        case S: return "S";
        case W: return "W";
        default: return "";
    }
}

char* format_dir(cardd dir) {
    char* c = malloc(3);

    strcat(c, get_dir(dir & N));
    strcat(c, get_dir(dir & S));
    strcat(c, get_dir(dir & W));
    strcat(c, get_dir(dir & E));

    return (c[0] == '\0') ? "0" : c;
}

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, cardd dir) {
    if (x < 0 || x > 2 || y < 0 || y > 2) {
        printf("Error: set_dir x=%d, y=%d - out of range\n", x, y);
        return;
    }
    else if (dir == (N|S) || dir == (W|E)) {
        printf("Error: set_dir dir=%s - direction not supported\n", format_dir(dir));
        return;
    }
    map[x][y] = dir;
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void) {
    for (int i=0; i<3; i++) {
        char c[3][3];
        for (int j=0; j<3; ++j) {
            strcpy(c[j], format_dir(map[i][j]));
        }
        printf("%-4s%s%4s\n", c[0], c[1], c[2]);
    }
}

int main(void) {
	// In dieser Funktion darf nichts verändert werden!
	set_dir(0, 1, N);
	set_dir(1, 0, W);
	set_dir(1, 4, W);
	set_dir(1, 2, E);
	set_dir(2, 1, S);

	show_map();

	set_dir(0, 0, N|W);
	set_dir(0, 2, N|E);
	set_dir(0, 2, N|S);
	set_dir(2, 0, S|W);
	set_dir(2, 2, S|E);
	set_dir(2, 2, E|W);

	show_map();

	return 0;
}
