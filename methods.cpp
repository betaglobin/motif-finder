#include "declaration.h"

char Full_Data::qual_file[] = "1";
char Full_Data::seq_file[] = "1";

void Full_Data::settings() {
    using namespace std;
    cout<<"==============================================="<<endl;
    cout<<"     Welcome to the Motif Searcher program"<<endl;
    cout<<"==============================================="<<endl;
    cout<<endl;
    cout<<"The following program finds if the sequences you provide contain any motifs"<<endl;
    cout<<"Motifs are created de novo"<<endl;
    cout<<endl;
    cout<<"STEP 1\n________________________________________________________\n";
    cout<<"Please choose the size of the motif you want to find..."<<endl;
    cout<<"Type a number between 4 and 7"<<endl;

    int check = 0;
    string answer;
    bool cond = true;
    std::string win_size;
    while (cond) {

        cin >> win_size;
        motif_len = std::stoi(win_size);

        if (motif_len < 4 || motif_len > 7) {
            cout << "Wrong motif size! Only numbers between 4 and 7 are allowed!" << endl;
            cout << "Type a number between 4 and 7:" << endl;
            check = 0;

        } else if (motif_len > 3 && motif_len < 8) {
            cout << "Window size selected: " << motif_len << endl;
            cout << "Do you want to keep it or change it?" << endl;
            cout << "Keep(1)/Change(2):" << endl;
            check = 1;
        }

        if(check == 1) {
            while (true) {
                cin >> answer;
                if (answer == "2") {
                    cout << "Ok." << endl;
                    cout << "Let's pick a new motif length:" << endl;
                    break;

                } else if (answer == "1") {
                    cout << "Confirmed..." << endl;
                    cond = false;
                    break;

                } else {
                    cout << "Choose 1 if you want to leave the length of the motif as it is" << endl;
                    cout << "Otherwise choose 2 to make a correction" << endl;
                }
            }
        }
    }
    reading_file();
};

void Full_Data::consider_qual(){
    using namespace std;
    cout<<"\nSTEP 3\n________________________________________________________\n";
    cout<<"Do you want to validate each nucleotide with the score?"<<endl;
    cout<<"Yes(1)/No(0)"<<endl;
    cout<<"By saying Yes YOU ARE AWARE that the lengths of the motifs can be reduce down to 4"<<endl;

    string answer;
    while (true){
        cin>>answer;

        if (answer == "1") {
            with_qual = 1;
            break;
        }

        else if (answer == "0") {
            with_qual = 0;
            break;
        }

        else
            cout<<"Wrong number... Choose 1 (with quality) or 0 (without)"<<endl;
    }

    cout<<"\nSTEP 4\n________________________________________________________\n";
    cout<<"Please provide a threshold of the confidence level"<<endl;
    cout<<"Values between 9-27"<<endl;

    while (true) {
        cin >> answer;

        int val = std::stoi(answer);

        if (val > 27) {
            cout << "Wrong number..." << endl;
        }
        else if (val < 9) {
            cout << "Wrong number..." << endl;
        }
        else {
            threshold = val;
            break;
        }
    }
};

void Full_Data::reading_file() {

    //the following methods reads content of the file into the one array of characters
    using namespace std;
    char file_name[20];
    bool cond = true;
    int position = 0;
    FILE *reading;
    int seq;

    while (cond) {
        cout<<"\nSTEP 2\n________________________________________________________\n";
        cout << "Provide a name of a sequence-containing file:" << endl;
        cin >> file_name;
        reading = fopen(file_name, "r");

        if (reading == NULL) {
            cout << "There is no such file!" << endl;

        }
        else {
            while ((seq = fgetc(reading)) != EOF) {
                seq_file[position] = seq;
                position += 1;
            }
            cond = false;
        }
        fclose(reading);
    }

    char qual_name[20];
    cond = true;
    position = 0;

    while (cond) {
        cout << "Provide a name of a quality scores-containing file:" << endl;
        cin >> file_name;
        reading = fopen(file_name, "r");

        if (reading == NULL) {
            cout << "There is no such file!" << endl;

        }
        else {
            while ((seq = fgetc(reading)) != EOF) {
                qual_file[position] = seq;
                position += 1;
            }
            cond = false;
        }
        fclose(reading);
    }
    consider_qual();
};

//==========================================================================================================
int Single_Seq::read_name(int begin) {

    //Single seq object is created dynamically and stored by Container object
    //it reads headers of the FASTA file, sequences and quality scores

    int i;
    int check = 0;
    bool cond = true;

    while (cond) {
        for (i = begin; i < FILESIZE; i++) { //start from the point you have finished last time
            char character = Full_Data::seq_file[i];

            if (character == '>') {
                start = i + 1;
            }
            else if (check == 0 && character == '\n'){ //if there is no blank character in the header take whole header as an ID
                end = i;
                stop = i;
                break;

            }
            else if (character == ' ' && check == 0) { //i'm interesting in first space only, and ID is everything what is on the left side of the space
                stop = i;
                check = 1;
            }
            else if (character == '\n') { //if it reaches end of the header
                end = i;
                break;
            }
        }

        std::string temp(&Full_Data::seq_file[start], &Full_Data::seq_file[stop]); //ID of the sequence
        name_string = temp;
        std::cout<<name_string<<"\t";
        cond = false;
    }
    return (end+1); //returns a positition from which another line will be examine
};

int Single_Seq::read_sequence(int begin){

    int i;
    end = begin;
    start = end;

    int check = 0;
    bool cond = true;

    while (cond) {
        for (i = start; i < FILESIZE; i++) { //start from the point you have finished last time
            char character = Full_Data::seq_file[i];

            if (character == '\n' && check == 0 &&Full_Data::seq_file[i+1] != '>') { //if a sequence has a '\n' inside of it
                gap = i;    //it DEALS ONLY with one or zero blank character within the sequence
                check = 1;
            }
            else if (character == '>') { //if another sequence's header appeared
                end = i - 1; //the sequence ends at that position
                break;
            }

            else if (character == '\0'){ //EOF has been reached
                end = i;
                break;
            }
        }

        if (check == 1) { //if there was a blank character take substrings...
            std::string part1(&Full_Data::seq_file[start], &Full_Data::seq_file[gap]);
            std::string part2(&Full_Data::seq_file[gap+1], &Full_Data::seq_file[end]);
            seq_string = part1 + part2; //...and concatenate them
            cond = false;
        }

        else if (check  == 0){
            std::string temp(&Full_Data::seq_file[start], &Full_Data::seq_file[end]);
            seq_string = temp;
            cond = false;
        }
    }

    std::cout << "Length of the sequence: " << seq_string.length() << std::endl;
    std::cout << seq_string << "\n" << std::endl;

    int size = int(seq_string.length());
    seq_array = new char[size];
    for (int i = 0; i < size; i++){
        char character = Full_Data::seq_file[start];
        start += 1;

        if (character == '\n'){
            i = i-1;
            continue;
        }
        seq_array[i] = character;
    }
    return end+1; //returns a position from which another line will be examine
};

int Single_Seq::read_quality(int begin, int length) {

    int i = 0;
    int check = 0;
    qual_vec.resize(length);

    for (int j = 0; j < length; j++) {

        char character2 = Full_Data::qual_file[begin+1];
        if (character2 == ' ' ||character2 == '\0'){
            check = 1;

        }
        else if(character2 != ' ' || character2 != '\0')
            check = 0;

        if (check == 1){
            std::string score(&Full_Data::qual_file[begin], &Full_Data::qual_file[begin + 1]);
            qual_vec[i] = score;
            begin +=2;
            i +=1;

        }
        else if (check == 0) {
            std::string score(&Full_Data::qual_file[begin], &Full_Data::qual_file[begin + 2]);
            qual_vec[i] = score;
            begin += 3;
            i += 1;
        }
    }
    return begin;
}

std::string Single_Seq::window_read(int begin, int end) {

    std::string temp_str = seq_string.substr(begin, end);
    return temp_str;
}

std::string Single_Seq::name_return(){

    return name_string;
}

std::string Single_Seq::sub_quality(int index) {

    return qual_vec[index];
}

//==========================================================================================================
void Container::filling_records(Full_Data object) {

    id_start = 0;
    qual_start = 0;
    jump = 0;
    int motif_len = object.motif_len;
    num_vertex = 0;

    lengths_of_seqs.resize(NUM_READS);

    for (int i = 0; i < NUM_READS; i++) {
        seq_start = record[i].read_name(id_start);

        qual_start = (seq_start - id_start) + jump;
        id_start = record[i].read_sequence(seq_start);

        int length = id_start-seq_start-2;
        num_vertex = num_vertex + length - motif_len+1;
        jump = record[i].read_quality(qual_start, length);
        lengths_of_seqs[i] = length;

    }
}

void Container::filling_vertex(Full_Data object) {
    int length = object.motif_len;
    int end = object.motif_len;
    int confidence = object.threshold;

    vertex_org.resize(num_vertex);
    int start = 0;
    filtered = 0;

    for (int i = 0; i < NUM_READS; i++){
        int begin = 0;
        std::string current_id = record[i].name_return();

        int iter = lengths_of_seqs[i] - length +1;

        for (int j = start; j < iter+start; j++){

            std::string temp_substr = record[i].window_read(begin, end);
            vertex_org[j].id_assigner(j, current_id, end);
            vertex_org[j].label_assigner(temp_substr);
            vertex_org[j].position_assigner(begin, length);

            int z = 0;
            for (int m = begin; m < length+begin; m++) {
                std::string score = record[i].sub_quality(m);
                vertex_org[j].score_assigner(z, score);
                z += 1;
            }

            int is_short = vertex_org[j].score_validation(length, confidence);
            vertex_org[j].print_vertex(length);
            begin += 1;

            if(is_short == 0){
                filtered += 1;

            }
        }
        start = iter+start;
    }

    vertex_mod.resize(filtered); //this is a number of qualified for being a motif-validated fragments
    int j = 0;
    for (int i = 0; i < filtered; i++){

        if (vertex_org[j].length_validation() == 1) {
            vertex_mod[i] = vertex_org[j];
            j+=1;
        }
        else if (vertex_org[j].length_validation() == 0){
            j+=1;
            i-=1;
        }
    }

    all_tags.resize(filtered+num_vertex);

    std::string header = vertex_org[0].return_name();
    std::string cmp_header = header;
    int last_ind = 0; int ind_o = 0; int ind_m = 0; int y = 0;

    while(ind_o < filtered+num_vertex){

        cmp_header = vertex_org[y].return_name();

        if (header == cmp_header) {
            all_tags[ind_o].id_seq = vertex_org[y].return_name();
            all_tags[ind_o].label = vertex_org[y].return_label(0);
            all_tags[ind_o].start = vertex_org[y].return_position(0);
            all_tags[ind_o].stop = vertex_org[y].return_position(1);
            ind_o+=1;
            last_ind = ind_o;
            y+=1;
        }

        else if (header != cmp_header){

            while(true) {
                cmp_header = vertex_mod[ind_m].return_name();

                if (header == cmp_header) {
                    all_tags[last_ind].id_seq = vertex_mod[ind_m].return_name();
                    all_tags[last_ind].label = vertex_mod[ind_m].return_label(1);
                    all_tags[last_ind].start = vertex_mod[ind_m].return_position(0);
                    all_tags[last_ind].stop = vertex_mod[ind_m].return_position(1);
                    ind_m+=1;
                    last_ind+=1;
                }

                else if (header != cmp_header){
                    ind_o = last_ind;
                    header = vertex_org[y].return_name();
                    break;
                }
            }
        }

    }

    for (int i = 0; i <filtered+num_vertex; i++){

        all_tags[i].id_num = i+1;
        all_tags[i].length = int(all_tags[i].label.length());
        std::cout<<all_tags[i].id_seq<<"    "<<all_tags[i].label<<"     "<<
                "len: "<<all_tags[i].length<<"    "<<all_tags[i].start<<"     "<<all_tags[i].stop<<"   "<<"id: "<<all_tags[i].id_num<<std::endl;
    }

}

int Vertex::return_position(int which) {
    if (which == 0){
        return pos_start;
    }
    else
        return pos_stop;
}

int Container::return_vertices(int which){

    if (which == 0)
        return num_vertex;
    else
        return filtered;
}

int Container::return_total_vs() {
    return num_vertex+filtered;
}

//==========================================================================================================
void Vertex::id_assigner(int uniq_num, std::string name, int score_len) {

    scores.resize(score_len);
    id_seq = name;
    id_number = uniq_num+1;
}

void Vertex::label_assigner(std::string substr) {

    label = substr;
}

void Vertex::position_assigner(int begin, int width) {

    pos_start = begin+1;
    pos_stop = pos_start+width;
}

void Vertex::score_assigner(int index, std::string score) {

    scores[index] = score;
}

void Vertex::print_vertex(int motif_len) {
    std::cout<<"Name: "<<id_seq<<"\t"<<label<<"\t"<<mod_label<<"\t"<<pos_start<<"\t"<<pos_stop<<"\tID: "<<id_number<< std::endl;
    for (int i = 0; i < motif_len; i++)
        std::cout<<scores[i]<<" ";
    std::cout<<"\n"<<std::endl;
}

int Vertex::length_validation(){

    if (mod_label.length() < 4)
        return 0;
    else
        return 1;
}

int Vertex::score_validation(int length, int threshold) {

    mod_label = label;
    std::queue<int> remove_pos;

    int j = 0;
    for (int i = 0; i < length; i++){
        std::string qual = scores[i];
        int int_val = std::stoi(qual);


        if (int_val < threshold){
            remove_pos.push(i);
            //std::cout<<"Score: "<<int_val<<'\t'<<"Position: "<<i<<std::endl;
            j+=1;
        }
    }

    int iter = 0;
    int curr_len = length;
    std::cout<<"\n"<<std::endl;
    while(!remove_pos.empty()){
        int index = ((remove_pos.front()) - iter);

        //std::cout<<mod_label<<" Current length: "<<curr_len<<" Position: "<<index<<std::endl;
        if(index == 0){
            mod_label = mod_label.substr(1, curr_len);
            iter +=1;
            curr_len-=1;
        }
        else if (index + 1 == curr_len-1){
            std::string last_el(1, mod_label[curr_len-1]);
            mod_label = mod_label.substr(0, index) + last_el;
            iter +=1;
            curr_len-=1;
        }

        else if (index != 0 && index != curr_len-1 && index+1 != curr_len-1){
            std::string part1 = mod_label.substr(0, index);
            std::string part2 = mod_label.substr(index+1, curr_len);
            mod_label = part1+part2;
            iter +=1;
            curr_len -=1;
        }

        else if (index == curr_len-1){
            mod_label = mod_label.substr(0, index);
            iter +=1;
            curr_len -=1;
        }

        remove_pos.pop();
    }

    if (mod_label.length() < 4){ //if the motif is shorter then 4, discard the whole object
        return 1;
    }
    else
        return 0;
}

std::string Vertex::return_name() {

    return id_seq;
}

int Vertex::return_id() {

    return id_number;
}

std::string Vertex::return_label(int which) {

    if (which == 0)
        return label;
    else
        return mod_label;

}
//==========================================================================================================

void Matrix::matrix_maker(Container container) {

    m_size = container.return_vertices(0) + 1;
    adj_matrix.resize(m_size);

    for (int i = 1; i < m_size; i++) {
        adj_matrix[i].resize(m_size);
    }

    for (int i = 1; i < m_size; i++) {
        for (int j = 1; j < m_size; j++)
            adj_matrix[i][j] = 0;
    }
}

void Matrix::connection(int v1, int v2){

    adj_matrix[v1][v2] = 1;
    adj_matrix[v2][v1] = 1;
}

void Matrix::edge_maker(Container container, Full_Data full_data) {

    if (full_data.with_qual == 1){

        int upper_bound = container.return_vertices(1);

        int id_nd, id_st;
        unsigned long len_st, len_nd;
        std::string seq_st, seq_nd, full_name_st, full_name_nd;

        for (int i = 0; i < upper_bound; i++){

            id_st = container.vertex_mod[i].return_id(); //id of the seq is the position in the graph of the first vertex in a pair
            seq_st = container.vertex_mod[i].return_label(1); //sequence to compare
            len_st = seq_st.length();
            full_name_st = container.vertex_mod[i].return_name();

            for (int j = i; j < upper_bound; j++) {
                seq_nd = container.vertex_mod[j].return_label(1);
                len_nd = seq_st.length();
                id_nd = container.vertex_mod[j].return_id();
                full_name_nd = container.vertex_mod[j].return_name();

                if (full_name_st != full_name_nd) { //examine only those that come from diffrent sequence

                    if (len_st == len_nd) {

                        if (seq_st == seq_nd) {
                            connection(id_st, id_nd);
                            std::cout << "Current ID is: " << id_st << "\tSeq: " << seq_st << " Len:\t" << len_st << " Name: "
                                      << full_name_st << std::endl;
                            std::cout << "Second ID is: " << id_nd << "\tSeq: " << seq_nd << " Len:\t" << len_nd << " Name: "
                                      << container.vertex_mod[j].return_name() << std::endl;
                        }
                    }
                    else if (len_st > len_nd) {

                        if (seq_st.find(seq_nd) != -1) {
                            id_nd = container.vertex_mod[j].return_id();
                            connection(id_st, id_nd);
                            std::cout << "Current ID is: " << id_st << "\tSeq: " << seq_st << " Len:\t" << len_st << " Name: "
                                      << full_name_st << std::endl;
                            std::cout << "Second ID is: " << id_nd << "\tSeq: " << seq_nd << " Len:\t" << len_nd << " Name: "
                                      << container.vertex_mod[j].return_name() << std::endl;
                        }
                    }
                    else if (len_st < len_nd) {

                        if (seq_nd.find(seq_st) != -1) {
                            id_nd = container.vertex_mod[j].return_id();
                            connection(id_st, id_nd);
                            std::cout << "Current ID is: " << id_st << "\tSeq: " << seq_st << " Len:\t" << len_st << " Name: "
                                      << full_name_st << std::endl;
                            std::cout << "Second ID is: " << id_nd << "\tSeq: " << seq_nd << " Len:\t" << len_nd << " Name: "
                                      << container.vertex_mod[j].return_name() << std::endl;
                        }
                    }
                }
            }
        }
    }

    else if (full_data.with_qual == 0) {

        int upper_bound = container.return_vertices(0);

        int id_nd, id_st;
        unsigned long len_st, len_nd;
        std::string seq_st, seq_nd, full_name;

        for (int i = 0; i < upper_bound; i++) {

            id_st = container.vertex_org[i].return_id(); //id of the seq is the position in the graph of the first vertex in a pair
            seq_st = container.vertex_org[i].return_label(0); //sequence to compare
            len_st = seq_st.length();
            full_name = container.vertex_org[i].return_name();

            std::cout << "Current ID is: " << id_st << "\tSeq: " << seq_st << " Len:\t" << len_st << " Name: "
                      << full_name << std::endl;

            for (int j = i; j < upper_bound; j++) {

                if (full_name != container.vertex_org[j].return_name()) { //examine only those that came from different sequence

                    seq_nd = container.vertex_org[j].return_label(0);
                    len_nd = seq_nd.length();
                    if (seq_st == seq_nd) {

                        id_nd = container.vertex_org[j].return_id();
                        connection(id_st, id_nd);
                        std::string full_name_nd = (container.vertex_org[j].return_name());

                        std::cout << "Second ID is: " << id_nd << "\tSeq: " << seq_nd << " Len:\t" << len_nd << " Name: "
                                  << full_name_nd << std::endl;
                    }
                }
            }
            std::cout<<"=============================================================="<<std::endl;
        }
    }
}

void Matrix::matrix_printer() {

    for (int i = 1; i < m_size; i++){
        for (int j = 1; j < m_size; j++){
            std::cout<<adj_matrix[i][j]<<"  ";
        }
        std::cout<<std::endl;
    }
}

void Matrix::clique_finder() {

    int size;
    int get_back = 0;
    int control = 0;

    for (int i = 1; i < m_size; i++) {
        std::vector<std::vector<int>> neigh_deep;
        std::vector<int> neigh_fixed;
        std::vector<int> subvec_size;
        std::queue<int> buffer;

        size = 1;

        for (int j = 1; j < m_size; j++) {

            if (adj_matrix[i][j] == 1) {             //load every adjacent vertex into the queue
                std::cout << i << "x" << j;
                std::cout << "Loading " << j << std::endl;
                buffer.push(j);
                neigh_fixed.resize(size);
                size += 1;
                control = 1;
            }
        }

        int z = 0;
        while (!buffer.empty()) {
            int v = buffer.front();
            neigh_fixed[z] = v;
            //std::cout << "Neigh fixed has: " << neigh_fixed[z] << std::endl;
            buffer.pop();
            z++;
        }

        if (control == 1) {

            size = size - 1;
            //std::cout << "SIZE IS: " << size << std::endl;
            //std::cout << "===============================" << std::endl;

            subvec_size.resize(size);
            neigh_deep.resize(size);

            for (int k = 0;
                 k < size; k++) {   //iterate size-times to initialize vectors in each cell of neigh_deep vector

                int dim_size = 1;
                int v = neigh_fixed[k]; //take the first neighbour and treat it as a zero-vertex
                neigh_deep[k].resize(dim_size);

                //std::cout << "INSPECTING : " << v << std::endl;

                for (int m = 1; m < m_size; m++) {

                    if (adj_matrix[v][m] == 1 && m != i) { //look for any adjacent vertices
                        //std::cout << "Vertex: " << v << " has: " << m << std::endl;
                        neigh_deep[k].resize(dim_size);
                        buffer.push(m);
                        dim_size += 1;

                    }
                }

                subvec_size[k] = dim_size - 1;
                //std::cout << "SUBVEC: " << k << " size: " << dim_size << std::endl;

                int index = 0;
                neigh_deep[k][index] = v; //zero element is a vertex "name" ==> length 1;

                //std::cout << neigh_deep[k][index] << " ";
                index += 1;

                while (!buffer.empty()) {
                    neigh_deep[k][index] = buffer.front();
                    //std::cout << neigh_deep[k][index] << " ";
                    buffer.pop();
                    index += 1;
                }

            }

            for (int l = 0; l < size; l++) {
                clique = 2;             //<=============== 9.01.2018
                if (get_back == 1)
                    break;

                if (subvec_size[l] > 1) { //i'm interesting only in inspecting neighbours that have ANY deep-neighbours
                    //std::cout << "I'm in cell l: " << neigh_fixed[l] << std::endl;
                    for (int n = 0; n < size; n++) {

                        if (get_back == 1)
                            break;

                        if (n != l) { //don't make any self-comparison of the vertex

                            if (subvec_size[n] > 1) { //and again inspect only those that have any deep-neighbours

                                std::cout << "Now i'm taking the n: " << n << std::endl;
                                for (int p = 1; p < subvec_size[l] + 1; p++) {

                                    if (neigh_deep[l][p] == neigh_deep[n][0]) {
                                        //std::cout << "Yes! I've found n: " << neigh_deep[n][0] << " in: p-cell number "
                                         //         <<
                                          //        p << " of l =" << neigh_fixed[l] << std::endl;

                                        clique += 1;
                                        //std::cout << "Clique is: " << clique << std::endl;
                                        clique_vertices.push(neigh_deep[n][0]);
                                        smaller_clique.push(neigh_deep[n][0]);
                                        //std::cout << (neigh_deep[n][0]) << " is loaded" << std::endl;
                                    }

                                    if (clique == 4) {
                                        clique_vertices.push(i);
                                        //std::cout << "Clique is: " << clique << std::endl;
                                        //std::cout << i << " is loaded" << std::endl;
                                        clique_vertices.push(neigh_deep[l][0]);
                                        //std::cout << (neigh_deep[l][0]) << " is loaded" << std::endl;
                                        get_back = 1;
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    if (get_back == 1) {
                        break;

                    }
                    else if (get_back == 0) {
                        clique = 2;
                        while (!clique_vertices.empty()){
                            clique_vertices.pop();
                        }
                    }
                }

                if (get_back == 1) {
                    break;
                }
            }
        }

        if (get_back == 1) {
            break;
        }

        if (control == 0) {
            while (!buffer.empty()) {
                buffer.pop();
            }
        }
    }
}

void Matrix:: motif_printer(Container container, Full_Data full_data) {
    using namespace std;

    int answer = full_data.with_qual;

    if (answer == 0) {

        int qv = clique_vertices.front() - 1;
        cout<<"\n\n====================================================================================="<<endl;
        cout << "Motif: " << container.vertex_org[qv].return_label(0)<<endl;

        while (!clique_vertices.empty()){
            qv = clique_vertices.front()-1;

            cout<<"seq: "<<container.vertex_org[qv].return_label(0)<<"\tID name: "<<container.vertex_org[qv].return_name()<<
                "\tstart: "<<container.vertex_org[qv].pos_start<<"\tstop: "<<
                container.vertex_org[qv].pos_stop<<"\tvertex num: "<< container.vertex_org[qv].return_id()<<endl;
            clique_vertices.pop();
        }
    }
    else if (answer == 1){

        std::vector<int> storage;
        std::queue<int> buffer;

        int c = 0;
        while(!clique_vertices.empty()){
            int v = clique_vertices.front();
            buffer.push(v);
            clique_vertices.pop();
            c += 1;
        }

        storage.resize(c);

        c = 0;

        while(!buffer.empty()){
            storage[c] = buffer.front();
            buffer.pop();
            c+=1;
        }

        string label = "m";
        string the_longest;
        unsigned long len = label.length();


        for (int i = 0; i < c; i++ ) {
            int qv = storage[i]-1;

            label = container.vertex_org[qv].return_label(1);
            //std::cout << label << endl;

            if (len < label.length()) {
                the_longest = label;
                len = the_longest.length();
            }
        }

        cout<<"\n\n====================================================================================="<<endl;
        cout <<"Full motif: "<<the_longest<<"\t All motifs: ";

        for (int j = 0; j < c; j++){
            int qv = storage[j]-1;
            cout << container.vertex_org[qv].return_label(1) << ", ";

        }
        cout<<endl;

        for (int i = 0; i < c; i++){
            int qv = storage[i]-1;

            cout<<"seq: "<<container.vertex_org[qv].return_label(1)<<"\tID name: "<<container.vertex_org[qv].return_name()<<
                "\t\tsomewhere_between: "<<container.vertex_org[qv].pos_start<<" - "<<container.vertex_org[qv].pos_stop<<endl;

        }
    }
}


//==========================================================================================================

void New_Matrix::matrix_maker(Container container) {

    m_size = container.return_total_vs()+1;
    adj_matrix.resize(m_size);

    for (int i = 1; i < m_size; i++) {
        adj_matrix[i].resize(m_size);
    }

    for (int i = 1; i < m_size; i++) {
        for (int j = 1; j < m_size; j++)
            adj_matrix[i][j] = 0;
    }

    for (int i = 1; i < m_size; i++) {
        for (int j = 1; j < m_size; j++)
            if(i == j)
                adj_matrix[i][j] = 1;
    }
}

void New_Matrix::connection(int v1, int v2){

    adj_matrix[v1][v2] = 1;
    adj_matrix[v2][v1] = 1;
}

void New_Matrix::edge_maker(Container container, Full_Data full_data) {


    int id_nd, id_st;
    unsigned long len_st, len_nd;
    std::string seq_st, seq_nd, full_name_st, full_name_nd;

    for (int i = 1; i < m_size-1; i++) {

        id_st = container.all_tags[i].id_num; //id of the seq is the position in the graph of the first vertex in a pair
        seq_st = container.all_tags[i].label; //sequence to compare
        len_st = seq_st.length();
        full_name_st = container.all_tags[i].id_seq;

        for (int j = i; j < m_size-1; j++) {
            id_nd = container.all_tags[j].id_num;
            seq_nd = container.all_tags[j].label;
            len_nd = seq_st.length();
            full_name_nd = container.all_tags[j].id_seq;

            if (full_name_st != full_name_nd) { //examine only those that come from different sequence

                if (len_st == len_nd) {

                    if (seq_st == seq_nd) {
                        connection(id_st, id_nd);
                        std::cout << "Current ID is: " << id_st << "\tSeq: " << seq_st << " Len:\t" << len_st
                                  << " Name: "
                                  << full_name_st << std::endl;
                        std::cout << "Second ID is: " << id_nd << "\tSeq: " << seq_nd << " Len:\t" << len_nd
                                  << " Name: "
                                  << full_name_nd << std::endl;
                    }
                } else if (len_st > len_nd) {

                    if (seq_st.find(seq_nd) != -1) {
                        id_nd = container.all_tags[j].id_num;
                        connection(id_st, id_nd);
                        std::cout << "Current ID is: " << id_st << "\tSeq: " << seq_st << " Len:\t" << len_st
                                  << " Name: "
                                  << full_name_st << std::endl;
                        std::cout << "Second ID is: " << id_nd << "\tSeq: " << seq_nd << " Len:\t" << len_nd
                                  << " Name: "
                                  << full_name_nd << std::endl;
                    }
                } else if (len_st < len_nd) {

                    if (seq_nd.find(seq_st) != -1) {
                        id_nd = container.all_tags[j].id_num;
                        connection(id_st, id_nd);
                        std::cout << "Current ID is: " << id_st << "\tSeq: " << seq_st << " Len:\t" << len_st
                                  << " Name: "
                                  << full_name_st << std::endl;
                        std::cout << "Second ID is: " << id_nd << "\tSeq: " << seq_nd << " Len:\t" << len_nd
                                  << " Name: "
                                  << full_name_nd << std::endl;
                    }
                }
            }
        }
    }
}






void New_Matrix::clique_finder() {

    using namespace std;

    multiplicate.resize(m_size);
    adj_mirror.resize(m_size);


    for (int i = 1; i < m_size; i++) {
        adj_mirror[i].resize(m_size);
        multiplicate[i].resize(m_size);
    }

    for (int j = 1; j < m_size; j++) {
        float sum = 0;

        for (int p = 1; p < m_size; p++) {
            if(adj_matrix[p][j] == 1)
                sum+=1;
        }
        if (sum != 0) {
            for (int i = 1; i < m_size; i++) {
                adj_matrix[i][j] = adj_matrix[i][j]/sum;
            }
        }
    }


    for (int i = 1; i < m_size; i++) {
        for (int j = 1; j < m_size; j++) {
            adj_mirror[i][j] = adj_matrix[i][j];

        }
    }

    //cout<<"PRZED PETLA"<<endl;

    for (int z = 0; z < 3; z++) {
        cout<<"PO PETLA"<<endl;

        for (int i = 1; i < m_size; i++)                // initializing elements...
            for (int j = 1; j < m_size; j++)
                multiplicate[i][j] = 0;


        for (int i = 1; i < m_size; i++) {
            for (int j = 1; j < m_size; j++) {
                for (int p = 1; p < m_size; p++) {
                    multiplicate[i][j] += adj_matrix[i][p] * adj_mirror[p][j];
                }
            }
        }

        for (int i = 1; i < m_size; i++) {
            for (int j = 1; j < m_size; j++) {
                adj_matrix[i][j] = multiplicate[i][j];
            }
        }

        for (int j = 1; j < m_size; j++) {
            for (int i = 1; i < m_size; i++) {
                adj_matrix[i][j] = (adj_matrix[i][j] * adj_matrix[i][j]);
            }
        }

        for (int j = 1; j < m_size; j++) {
            for (int i = 1; i < m_size; i++) {
                if (adj_matrix[i][j] < 0.01)
                    adj_matrix[i][j] = 0;
                if (adj_matrix[i][j] > 0.99) {
                    adj_matrix[i][j] = 1;

                }
            }
        }

        for (int j = 1; j < m_size; j++) {
            float sum = 0;

            for (int p = 1; p < m_size; p++) {
                sum = sum + adj_matrix[p][j];
            }
            for (int i = 1; i < m_size; i++) {
                if (sum != 0)
                    adj_matrix[i][j] = adj_matrix[i][j] / sum;
            }
        }

        for (int m = 1; m < m_size; m++) {
            for (int n = 1; n < m_size; n++) {
                adj_mirror[m][n] = adj_matrix[m][n];
            }
        }
    }
}

void New_Matrix:: motif_printer(Container container, Full_Data full_data) {
    using namespace std;

    cout<<"MOTYW"<<endl;
    int the_largest = 0;
    int current = 0;
    int where = 1;
    for (int m = 1; m < m_size; m++) {
        for (int n = 1; n < m_size; n++) {
            if (adj_matrix[m][n] > 0) {
                current += 1;
                //cout<<"JEST JEDYNKA"<<endl;
            }

        }

        if (current > the_largest) {
            the_largest = current;
            where = m;
        }

        current = 0;
    }

    queue<int> buffer;

    for (int j = 1; j < m_size; j++) {
        if (adj_matrix[where][j] > 0) {
            buffer.push(j);
        }
    }

    while (!buffer.empty()) {
        int qv = buffer.front() - 1;
        cout << "seq: " << container.all_tags[qv].label << "\tID name: " << container.all_tags[qv].id_seq <<
             "\tstart: " << container.all_tags[qv].start << "\tstop: " <<
             container.all_tags[qv].stop << "\tvertex num: " << container.all_tags[qv].id_num << endl;
        buffer.pop();
    }
}
