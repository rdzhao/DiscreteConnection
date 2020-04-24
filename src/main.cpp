#include "Connection.h"

int main(int argc, char* argv[]){
    Connection con;
    con.read(argv[1]);
    con.preprocessing();
    con.compute();
    con.write();

    return 1;
}