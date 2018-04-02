#include "declaration.h"

int main() {
    Full_Data Reader;
    Reader.settings();

    Container Shelf;
    Shelf.initializer();
    Shelf.filling_records(Reader);
    Shelf.filling_vertex(Reader);

    Matrix Graph;
    Graph.matrix_maker(Shelf);
    Graph.edge_maker(Shelf, Reader);
    Graph.clique_finder();
    Graph.motif_printer(Shelf, Reader);



    New_Matrix New_Graph;
    New_Graph.matrix_maker(Shelf);
    New_Graph.edge_maker(Shelf, Reader);
    New_Graph.clique_finder();
    New_Graph.motif_printer(Shelf, Reader);

    return 0;
}