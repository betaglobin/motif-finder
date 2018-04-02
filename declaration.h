#ifndef METHODS_H_
#define METHODS_H_

#include <iostream>
#include <string>
#include <fstream>
#include <limits>
#include <new>
#include <vector>
#include <queue>

#define FILESIZE 20000
#define NUM_READS 7


class Full_Data{
    friend class Container;
    friend class Matrix;

public:
    void settings();
    void reading_file();
    void consider_qual();

    static char qual_file[FILESIZE];
    static char seq_file[FILESIZE];


private:
    int motif_len;
    int threshold;
    int with_qual;
};


class Single_Seq : public Full_Data{
    friend class Full_Data;

public:
    int read_sequence(int begin);
    int read_name(int begin);
    int read_quality(int begin, int length);
    std::string window_read(int start, int stop);
    std::string name_return();
    std::string sub_quality(int index);

private:
    std::vector <std::string> qual_vec;

    std::string seq_string;

    std::string name_string;

    char* name_array;
    char* seq_array;

    int start;
    int stop;
    int end;
    int gap;

};

class Vertex{
    friend class Matrix;
    friend class Container;

public:
    void id_assigner(int uniq_num, std::string name, int score_len);
    void label_assigner(std::string substr);
    void position_assigner(int begin, int width);
    void score_assigner(int index, std::string score);

     int score_validation(int length, int threshold);
     int length_validation();

     std::string return_name();
     std::string return_label(int which);
     int return_id();
     int return_position(int which);

    void print_vertex(int motif_len);

private:
    std::string id_seq;
    std::string label;
    std::string mod_label;
    int pos_start;
    int pos_stop;
    int id_number;

    std::vector <std::string> scores;
};

class Tag{
public:
    std::string label;
    std::string id_seq;
    int start;
    int stop;
    int id_num;
    int length;


};


class Container : public Full_Data{
    friend class Single_Seq;
    friend class Full_Data;

public:
    std::vector<Single_Seq> record;
    std::vector<Vertex> vertex_org;
    std::vector<Vertex> vertex_mod;

    std::vector<int> lengths_of_seqs;

    std::vector <Tag> all_tags;

    void initializer(){record.resize(NUM_READS);}
    void filling_records(Full_Data object);
    void filling_vertex(Full_Data object);

    int return_vertices(int which);

    int return_total_vs();


private:
    int id_start;
    int seq_start;
    int qual_start;
    int jump;
    int num_vertex;     // number of unmodified vertices
    int filtered;       //number of score-validated vertices

};


class Matrix{
    friend class Container;
    friend class Full_Data;
    friend class Vertex;

public:
    void matrix_maker(Container container);
    void edge_maker(Container container, Full_Data full_data);
    void connection(int v1, int v2);
    void clique_finder();
    void matrix_printer();
    void motif_printer(Container container, Full_Data full_data);


private:
    std::vector<std::vector <int> > adj_matrix;
    int m_size;

    int clique;
    int max_clique;
    std::queue<int>clique_vertices;
    std::queue<int>smaller_clique;

};

class New_Matrix{
    friend class Container;
    friend class Full_Data;
    friend class Vertex;

public:
    void matrix_maker(Container container);
    void edge_maker(Container container, Full_Data full_data);
    void connection(int v1, int v2);
    void clique_finder();
    void matrix_printer();
    void motif_printer(Container container, Full_Data full_data);


private:
    std::vector<std::vector <float> > adj_matrix;
    std::vector<std::vector <float> > multiplicate;
    std::vector<std::vector <float> > adj_mirror;
    int m_size;

    int clique;
    int max_clique;
    std::queue<int>clique_vertices;
    std::queue<int>smaller_clique;

};

#endif
