
#include "bclist_parallel.h"
#include <boost/math/special_functions/math_fwd.hpp>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <omp.h>
#include <sys/types.h>
#include <utility>
#include <vector>
#include <boost/math/special_functions/binomial.hpp>

double time_counter0 = 0.0;
double time_counter1 = 0.0;
lint combination_cache[MAX_N+1][MAX_K+1] = {0};
lint my_factorial_cache[MAX_N+1][MAX_K+1] = {0};

VertexDegree::VertexDegree(): vertex(-1), degree(-1){}

VertexDegree::VertexDegree(unsigned v, unsigned d): vertex(v), degree(d){}

VertexDegree::~VertexDegree(){}

vector<unsigned> Tools::intersection(vector<unsigned> &vec1, vector<unsigned> &vec2){
    return intersection(vec1, vec2, 0, 0);
}

vector<unsigned> Tools::intersection(vector<unsigned> &vec1, vector<unsigned> &vec2, int offset1, int offset2){
    vector<unsigned> ans;
    if(vec1.empty() || vec2.empty() || offset1 >= vec1.size() || offset2 >= vec2.size()){
        return ans;
    }
    
    // we assume vec1 and vec2 are sorted by ascending order of the elements
    unsigned i, j;
    i = offset1;
    j = offset2;
    while(i < vec1.size() && j < vec2.size()){
        if(vec1[i] == vec2[j]){
            ans.push_back(vec1[i]);
            i++;
            j++;
        }
        else if(vec1[i] < vec2[j]){
            i++;
        }
        else{
            j++;
        }
    }
    return ans;
}

unsigned Tools::intersection_count(vector<unsigned> &vec1, vector<unsigned> &vec2){
    unsigned ans = 0;
    if(vec1.empty() || vec2.empty()){
        return 0;
    }
    
    // we assume vec1 and vec2 are sorted by ascending order of the elements
    unsigned i, j;
    i = j = 0;
    while(i < vec1.size() && j < vec2.size()){
        if(vec1[i] == vec2[j]){
            ans++;
            i++;
            j++;
        }
        else if(vec1[i] > vec2[j]){
            j++;
        }
        else{
            i++;
        }
    }
    
    return ans;
}

lint Tools::choose(lint n, lint k){
    if(combination_cache[n][k] > 0){
        return combination_cache[n][k];
    }
    else{
        if (k > n) {
            return 0;
        }
        lint backup_n = n;
        lint backup_k = k;
        lint r = 1;
        for (lint d = 1; d <= k; ++d) {
            r *= n--;
            r /= d;
        }
        if(n <= MAX_N && k <= MAX_K){
            combination_cache[backup_n][backup_k] = r;
        }
        return r;
    }
}

lint Tools::my_factorial(lint n, lint k){
    if(my_factorial_cache[n][k] > 0){
        return my_factorial_cache[n][k];
    }
    else{
        if(k > n){
            return 0;
        }
        lint r = 1;
        for(lint i = k; i <= n; ++i){
            r *= i;
        }
        if(n <= MAX_N && k <= MAX_K){
            my_factorial_cache[n][k] = r;
        }
        return r;
    }
}

// below are special bigraphs

SpecialBigraph::SpecialBigraph(string path, unsigned p_value, unsigned q_value, unsigned priority_strategy){
    clear_everything();
    data_file_path = path;
    p = p_value;
    q = q_value;
    priority = priority_strategy;
}

void SpecialBigraph::clear_everything(){
    largest_index_in_partition[0] = largest_index_in_partition[1] = 0;
    num_vertices = 0;
    num_edges = 0;
    result_count = 0;
    bcliques.clear();
    edges.clear();
    vertices[0].clear(); vertices[1].clear();
    vertices_in_left.clear();
    vertices_in_right.clear();
    adj_vec.clear();
}

bool SpecialBigraph::all_num(string &s) {
    for (int i = 0; i < SZ(s); i++) if ((s[i] >= '0' && s [i] <= '9') == false) return false;
    return true;
}

void SpecialBigraph::add_vertex(unsigned A, unsigned side){
    if(vertices[side].find(A) == vertices[side].end()){
        if(side == 0)
            vertices_in_left.emplace_back(A);
        else
            vertices_in_right.emplace_back(A);
        vertices[side][A] = 1;
    }
}

void SpecialBigraph::add_edge(unsigned &A, unsigned &B){
    add_vertex(A, 0);
    add_vertex(B, 1);
    if(edges.find(make_pair(A, B)) == edges.end()){
        edges.insert(make_pair(A, B));
        num_edges++;
    }
}

void SpecialBigraph::get_index(unsigned &A, unsigned side){
    if(vertices[side].find(A) == vertices[side].end()){
        vertices[side][A] = largest_index_in_partition[side]++;
    }
    A = vertices[side][A];
}

void SpecialBigraph::read_graph(){

    std::ifstream aaa;
    aaa.open(data_file_path.c_str(),std::ios::in);

    string s;
    while(getline(aaa, s)){
        stringstream ss;
        ss << s;
        vector <string> vec_str;
        for(string z; ss >> z; vec_str.push_back(z));
        if(SZ(vec_str) >= 2){
            bool is_all_number = true;
            for(unsigned i = 0; i < MIN(2, SZ(vec_str)); i++)
                is_all_number &= all_num(vec_str[i]);
            if(is_all_number){
                unsigned A, B;
                ss.clear(); ss << vec_str[0]; ss >> A;
                ss.clear(); ss << vec_str[1]; ss >> B;
                add_edge(A, B);
            }
        }
    }
    vertices[0].clear();
    vertices[1].clear();
    
    n_vertices[0] = vertices_in_left.size();
    n_vertices[1] = vertices_in_right.size();
    num_vertices = n_vertices[0] + n_vertices[1];
    num_edges = edges.size();
    deg = new unsigned[num_vertices]();
    adj_vec.resize(num_vertices, vector<unsigned>());
    
    largest_index_in_partition[0] = 0;
    largest_index_in_partition[1] = n_vertices[0];
    for(set < pair <unsigned, unsigned> >::iterator it = edges.begin(); it != edges.end(); it++){
        unsigned A = it->first;
        unsigned B = it->second;
        get_index(A, 0);
        get_index(B, 1);
        adj_vec[A].emplace_back(B);
        adj_vec[B].emplace_back(A);
        list_of_edges.emplace_back(make_pair(A, B));
        deg[A]++;
        deg[B]++;
    }
    
    #ifdef IS_DEBUGGING
    //print_adj();
    
    //print_edges();
    #endif
    #ifndef FREESUB
    printf("end of read graph\n");
    printf("#vertices = %u, #left_vertices = %u, #right_vertices = %u",num_vertices,n_vertices[0],n_vertices[1]);
    printf("#edges = %d\n",num_edges);
    #endif
}



void SpecialBigraph::trim_graph_by_core(){
    //cout << p << " " << q << endl;
    
    unsigned start_idx, end_idx;
    start_idx = end_idx = 0;
    
    unsigned *to_remove_vertices = new unsigned[num_vertices];
    bool *removed = new bool[num_vertices];
    
    for(unsigned i = 0; i < num_vertices; i++){
        to_remove_vertices[i] = 0;
        removed[i] = false;
    }
    
    // step 1: collect the initial vertices with degree less than q in left and p in right
    for(unsigned i = 0; i < num_vertices; i++){
        if((i < n_vertices[0] && deg[i] < q) || (i >= n_vertices[0] && deg[i] < p)){
            to_remove_vertices[end_idx++] = i;
            removed[i] = true;
        }
    }
        
    // step 2: recursively remove all vertices with degree less than q in left and p in right, i.e., (q,p)-core
    while(start_idx != end_idx){
        unsigned vertex = to_remove_vertices[start_idx++];
        
        //cout << "remove : " << vertex << endl;
        for(int i = 0; i < adj_vec[vertex].size(); i++){
            unsigned other = adj_vec[vertex][i];
            if(!removed[other]){
                deg[other]--;
                
                if((other < n_vertices[0] && deg[other] < q) || (other >= n_vertices[0] && deg[other] < p)){
                    to_remove_vertices[end_idx++] = other;
                    removed[other] = true;
                }
            }
        }
        adj_vec[vertex].clear();
        deg[vertex] = 0;
    }
        
    for(unsigned i = 0; i < num_vertices; i++){
        if(!removed[i]){
            vector<unsigned> new_vec(deg[i]);
            unsigned idx = 0;
            for(unsigned j = 0; j < adj_vec[i].size(); j++){
                if(!removed[adj_vec[i][j]]){
                    new_vec[idx++] = adj_vec[i][j];
                }
            }
            adj_vec[i].clear();
            adj_vec[i] = new_vec;
        }
    }
    
    #ifdef IS_DEBUGGING
    print_adj();
    #endif
    
    // printf("after trimming by core\n");
    
    reformat_graph();
    
    // printf("end of trim graph by core\n");
}

void SpecialBigraph::trim_graph_by_two_hop(){
    
    unsigned array_size = n_vertices[0];
    
    unsigned *common_neig_count = new unsigned[array_size]();
    unsigned *common_neig_map = new unsigned[array_size]();
    unsigned *aux_array_two_neig = new unsigned[array_size]();
    for(unsigned i = 0; i < array_size; i++){
        unsigned idx = 0;
        for(unsigned j = 0; j < adj_vec[i].size(); j++){
            unsigned neighbor = adj_vec[i][j];
            for(unsigned k = 0; k < adj_vec[neighbor].size(); k++){
                unsigned two_hop_neighbor = adj_vec[neighbor][k];
                if(two_hop_neighbor < i){
                    common_neig_map[two_hop_neighbor]++;
                    if(common_neig_map[two_hop_neighbor] == 1){
                        aux_array_two_neig[idx++] = two_hop_neighbor;
                    }
                }
                else{
                    break;
                }
            }
        }
        for(unsigned j = 0; j < idx; j++){
            unsigned two_hop_neighbor = aux_array_two_neig[j];
            if((i < n_vertices[0] && common_neig_map[two_hop_neighbor] >= q) || (i >= n_vertices[0] && common_neig_map[two_hop_neighbor] >= p)){
                common_neig_count[i]++;
                common_neig_count[two_hop_neighbor]++;
            }
            common_neig_map[two_hop_neighbor] = 0;
        }
    }
    
    printf("count common neighbors\n");
    
    for(unsigned i = 0; i < array_size; i++){
        if((i < n_vertices[0] && common_neig_count[i] < p - 1) || (i >= n_vertices[0] && common_neig_count[i] < q - 1)){
            for(unsigned j = 0; j < deg[i]; j++){
                unsigned other = adj_vec[i][j];
                for(unsigned k = 0; k < deg[other]; k++){
                    if(adj_vec[other][k] == i){
                        adj_vec[other][k] = adj_vec[other][--deg[other]];
                        break;
                    }
                }
            }
            adj_vec[i].clear();
            deg[i] = 0;
        }
    }
    delete[] common_neig_count;
    delete[] common_neig_map;
    delete[] aux_array_two_neig;
    
    #ifdef IS_DEBUGGING
    print_adj();
    #endif
    
    reformat_graph();
    
    printf("end of trim graph by two hop neighbors\n");
}

void SpecialBigraph::reformat_graph(){
    // re-format the bipartite graph since we might have removed some vertices
    edges.clear();
    for(unsigned i = 0; i < n_vertices[0]; i++){
        for(unsigned j = 0; j < deg[i]; j++){
            edges.insert(make_pair(i, adj_vec[i][j]));
        }
    }
    num_edges = edges.size();
    
    unsigned removed_count[2];
    removed_count[0] = removed_count[1] = 0;
    unsigned side;
    for(unsigned i = 0; i < num_vertices; i++){
        if(deg[i] == 0){
            side = i < n_vertices[0] ? 0 : 1;
            removed_count[side]++;
        }
    }
    
    n_vertices[0] = n_vertices[0] - removed_count[0];
    n_vertices[1] = n_vertices[1] - removed_count[1];
    
    num_vertices = n_vertices[0] + n_vertices[1];
    
    list_of_edges.clear();
    for(unsigned i = 0; i < adj_vec.size(); i++){
        adj_vec[i].clear();
    }
    adj_vec.clear();
    
    //cout << num_vertices << endl;
    
    vertices[0].clear();
    vertices[1].clear();
    delete[] deg;
    deg = new unsigned[num_vertices];
    memset(deg, 0, num_vertices * sizeof(unsigned));
    adj_vec.resize(num_vertices, vector<unsigned>());
    
    largest_index_in_partition[0] = 0;
    largest_index_in_partition[1] = n_vertices[0];
    for(set < pair <unsigned, unsigned> >::iterator it = edges.begin(); it != edges.end(); it++){
        
        unsigned A = it->first;
        unsigned B = it->second;
        
        //cout << A << " " << B << endl;
        
        get_index(A, 0);
        get_index(B, 1);
        
        //cout << "new : " << A << " " << B << endl;
        
        adj_vec[A].emplace_back(B);
        adj_vec[B].emplace_back(A);
        list_of_edges.emplace_back(make_pair(A, B));
        deg[A]++;
        deg[B]++;
    }
    
    for(unsigned i = 0; i < num_vertices; i++){
        sort(adj_vec[i].begin(), adj_vec[i].end());
    }
    
    #ifdef IS_DEBUGGING
    print_adj();
    
    print_edges();
    #endif
    #ifndef FREESUB
    printf("end of reformat graph\n");
    printf("#vertices = %u, #left_vertices=%u, #right_vertices=%u\n",num_vertices,n_vertices[0],n_vertices[1]);
    printf("#edges = %lu\n",edges.size());
    #endif  
}

void SpecialBigraph::collect_two_hop_adj(){
    unsigned array_size = n_vertices[0];
    
    two_hop_adj_maxsize = new unsigned[array_size];
    two_hop_adj_vec = new unsigned*[array_size];
    for(unsigned i = 0; i < array_size; i++){
        two_hop_adj_maxsize[i] = NTWOHOPS;
        two_hop_adj_vec[i] = new unsigned[two_hop_adj_maxsize[i]];
    }
    two_hop_adj_size = new unsigned[array_size]();
    
    unsigned *common_neig_map = new unsigned[array_size]();
    unsigned *aux_array_two_neig = new unsigned[array_size]();
    
    for(unsigned i = 0; i < array_size; i++){
        unsigned idx = 0;
        for(unsigned j = 0; j < adj_vec[i].size(); j++){
            unsigned neighbor = adj_vec[i][j];
            for(unsigned k = 0; k < adj_vec[neighbor].size(); k++){
                unsigned two_hop_neighbor = adj_vec[neighbor][k];
                if(two_hop_neighbor < i){
                    common_neig_map[two_hop_neighbor]++;
                    if(common_neig_map[two_hop_neighbor] == 1){
                        aux_array_two_neig[idx++] = two_hop_neighbor;
                    }
                }
                else{
                    break;
                }
            }
        }
        
        for(unsigned j = 0; j < idx; j++){
            unsigned two_hop_neighbor = aux_array_two_neig[j];
            if(common_neig_map[two_hop_neighbor] >= q){
                if(two_hop_adj_size[two_hop_neighbor] >= two_hop_adj_maxsize[two_hop_neighbor]){
                    //two_hop_adj_maxsize[two_hop_neighbor] += NTWOHOPS;
                    two_hop_adj_maxsize[two_hop_neighbor] *= 2;
                    unsigned *temp_array = new unsigned[two_hop_adj_maxsize[two_hop_neighbor]];
                    memcpy(temp_array, two_hop_adj_vec[two_hop_neighbor], sizeof(unsigned)*two_hop_adj_size[two_hop_neighbor]);
                    delete[] two_hop_adj_vec[two_hop_neighbor];
                    two_hop_adj_vec[two_hop_neighbor] = temp_array;
                }
                two_hop_adj_vec[two_hop_neighbor][two_hop_adj_size[two_hop_neighbor]++] = i;
            }
            common_neig_map[two_hop_neighbor] = 0;
        }
    }

    delete[] common_neig_map;
    delete[] aux_array_two_neig;
}
    
bool SpecialBigraph::prepare_graph(){
    //标记被prepare了，则准备析构
    if_prepare=1;

    #ifdef USE_CORE_REDUCTION
    trim_graph_by_core();
    #endif
    
    //trim_graph_by_two_hop();
    //print_deg();
    
    if(n_vertices[0] == 0 || n_vertices[1] == 0){
        #ifndef FREESUB 
        printf("No results because the graph is pruned by core\n");
        #endif
        return false;
    }
    
    if(estimate_cost(0) < estimate_cost(1)){
        anchor_left = true;
    }
    else{
        anchor_left = false;
    }


    #ifndef FREESUB 
    printf("anchor_left = %u\n",anchor_left);
    #endif
    
    sort_vertices(priority);
    
    #ifdef IS_DEBUGGING
    print(all_vertices);
    #endif

    #ifndef FREESUB 
    printf("finish sorting vertices\n");
    #endif
    
    unsigned id = 0;
    for(vector<unsigned>::iterator it = all_vertices.begin(); it != all_vertices.end(); it++){
        vertexids[*it] = id++;
    }
    
    // sort the graph by the rank of vertices
    for(unsigned i = 0; i < adj_vec.size(); i++){
        adj_vec[i].clear();
    }
    adj_vec.clear();
    delete[] deg;
    deg = new unsigned[num_vertices]();
    adj_vec.resize(num_vertices, vector<unsigned>());
    for(vector< pair<unsigned, unsigned> >::iterator it = list_of_edges.begin(); it != list_of_edges.end(); it++){
        unsigned A = vertexids[it->first];
        unsigned B = vertexids[it->second];

        adj_vec[A].emplace_back(B);
        adj_vec[B].emplace_back(A);
        //list_of_edges.push_back(make_pair(A, B));
        deg[A]++;
        deg[B]++;
    }
    
    if(!anchor_left){
        unsigned tmp_num = n_vertices[0];
        n_vertices[0] = n_vertices[1];
        n_vertices[1] = tmp_num;
        
        unsigned tmp_value = p;
        p = q;
        q = tmp_value;
    }
    
    for(unsigned i = 0; i < num_vertices; i++){
        sort(adj_vec[i].begin(), adj_vec[i].end());
    }
    
    collect_two_hop_adj();
    
    for(unsigned i = 0; i < n_vertices[0]; i++){
        sort(two_hop_adj_vec[i], two_hop_adj_vec[i] + two_hop_adj_size[i]);
    }
    
    #ifdef IS_DEBUGGING
    cout << "adj after preparing" << endl;
    print_adj();
    
    print_two_hop_adj();
    #endif
    
    return true;
}


double SpecialBigraph::estimate_cost(unsigned side){
    srand (time(NULL));
    
    unsigned num_rounds = num_vertices * SAMPLE_RATE;
    
    lint total_two_hop_deg = 0;
    lint max_two_hop_deg = 0;
    
    unsigned *common_neig_map = new unsigned[num_vertices]();
    unsigned *aux_array_two_neig = new unsigned[num_vertices]();
    unsigned offset = side == 0 ? 0 : n_vertices[0];
    unsigned common_neig_threshold = side == 0 ? q : p;
    
    for(unsigned r = 0; r < num_rounds; r++){
        lint estimated_two_hop_deg = 0;
        unsigned i = rand() % n_vertices[side] + offset;
        unsigned idx = 0;
        for(unsigned j = 0; j < adj_vec[i].size(); j++){
            unsigned neighbor = adj_vec[i][j];
            for(unsigned k = 0; k < adj_vec[neighbor].size(); k++){
                unsigned two_hop_neighbor = adj_vec[neighbor][k];
                if(two_hop_neighbor < i){
                    common_neig_map[two_hop_neighbor]++;
                    if(common_neig_map[two_hop_neighbor] == 1){
                        aux_array_two_neig[idx++] = two_hop_neighbor;
                    }
                }
                else{
                    break;
                }
            }
        }
    
        for(unsigned j = 0; j < idx; j++){
            unsigned two_hop_neighbor = aux_array_two_neig[j];
            if(common_neig_map[two_hop_neighbor] >= common_neig_threshold){
                estimated_two_hop_deg += 1;
            }
            common_neig_map[two_hop_neighbor] = 0;
        }
        
        max_two_hop_deg = max_two_hop_deg < estimated_two_hop_deg ? estimated_two_hop_deg : max_two_hop_deg;
        total_two_hop_deg = total_two_hop_deg + estimated_two_hop_deg * n_vertices[side];
    }
    
    try {
        if(num_rounds == 0 ) throw "123" ;
        total_two_hop_deg = total_two_hop_deg / num_rounds;
    } catch (...) {
       return 0;
    }

    // total_two_hop_deg = total_two_hop_deg / num_rounds;
    #ifndef FREESUB
    cout << "estimated_total_two_hop_deg[" << side << "]=" << total_two_hop_deg << ", estimated_avg_two_hop_deg[" << side << "]=" << total_two_hop_deg / n_vertices[side] << ", max_two_hop_deg[" << side << "]=" << max_two_hop_deg << endl;
    #endif  
    
    lint avg_two_hop_deg = MAX(2, total_two_hop_deg / n_vertices[side]);  // we let avg_two_hop_deg be at least 2 to deal with corner case
    lint pq_value = side == 0? p:q;
    double totalCost = total_two_hop_deg * pow(avg_two_hop_deg, pq_value-2);

    #ifndef FREESUB
    cout << "totalCost = " << totalCost << endl;
    #endif  
    
    delete[] common_neig_map;
    delete[] aux_array_two_neig;
    return totalCost;
}

void SpecialBigraph::find_combinations(vector< vector<unsigned> > &combs, lint &comb_count, vector<unsigned> &seed_vertices, vector<unsigned> &comb, unsigned beg_offset, unsigned curr_depth, unsigned total_depth){
    //comb用于临时存储结果。len(comb)==total_depth；beg_offset为左侧游标，初始值取0；
    // total_depth是取出个数；curr_depth用于指示递归深度，初始值取total_depth）
    unsigned N = seed_vertices.size();
    if (curr_depth == 0) {
        #ifndef COUNT_ONLY
        combs.push_back(comb);
        #endif
        comb_count++;
        return;
    }
    for (unsigned i = beg_offset; i < N; i++){
        comb[total_depth-curr_depth] = seed_vertices[i];
        find_combinations(combs, comb_count, seed_vertices, comb, i + 1, curr_depth - 1, total_depth);
    }
}

void SpecialBigraph::sort_vertices(unsigned strategy){
    switch (strategy) {
        case 0: {// random vertex order
            // printf("random vertex order is applied...\n"); ;
            
            all_vertices.resize(num_vertices);
            for(unsigned i = 0; i < num_vertices; i++){
                all_vertices[i] = i;
            }
            
            break;
        }
        case 1: {// degree vertex order
            // printf("degree vertex order is applied...\n"); ;
            
            //cout << "num_vertices = " << num_vertices << endl;
            
            vector<VertexDegree> verdegs(num_vertices);
            for(unsigned i = 0; i < num_vertices; i++){
                VertexDegree vd(i, deg[i]);
                verdegs[i] = vd;
            }
            sort(verdegs.begin(), verdegs.end());
            
            all_vertices.resize(num_vertices);
            for(unsigned i = 0; i < num_vertices; i++){
                all_vertices[i] = verdegs[i].vertex;
            }
            
            break;
        }
        case 2:{ // core vertex order:
            // printf("core vertex order is applied...\n"); ;
            
            unsigned *to_remove_vertices = new unsigned[num_vertices];
            bool *removed = new bool[num_vertices];
            
            all_vertices.resize(num_vertices);
            unsigned idx = num_vertices;
            for(unsigned i = 0; i < num_vertices; i++){
                to_remove_vertices[i] = 0;
                removed[i] = false;
            }
            
            unsigned p_value = p;
            unsigned q_value = q;
            
            unsigned start_idx, end_idx;
            start_idx = end_idx = 0;
            while(idx > 0){
                for(unsigned i = 0; i < end_idx; i++){
                    to_remove_vertices[i] = 0;
                }
                start_idx = end_idx = 0;
                
                for(unsigned i = 0; i < num_vertices; i++){
                    if(((i < n_vertices[0] && deg[i] < q_value) || (i >= n_vertices[0] && deg[i] < p_value)) && !removed[i]){
                        to_remove_vertices[end_idx++] = i;
                        removed[i] = true;
                        all_vertices[--idx] = i;
                    }
                }
                
                while(start_idx != end_idx){
                    unsigned vertex = to_remove_vertices[start_idx++];
     
                    for(int i = 0; i < deg[vertex]; i++){
                        unsigned other = adj_vec[vertex][i];
                        
                        for(int j = 0; j < deg[other]; j++){
                            if(adj_vec[other][j] == vertex){
                                unsigned tmp = adj_vec[other][j];
                                adj_vec[other][j] = adj_vec[other][--deg[other]];
                                adj_vec[other][deg[other]] = tmp;
                                break;
                            }
                        }
                        
                        if(((other < n_vertices[0] && deg[other] < q_value) || (other >= n_vertices[0] && deg[other] < p_value)) && !removed[other]){
                            to_remove_vertices[end_idx++] = other;
                            removed[other] = true;
                            all_vertices[--idx] = other;
                        }
                    }
                    deg[vertex] = 0;
                }
                p_value++;
                q_value++;
                
            }
            delete[] to_remove_vertices;
            delete[] removed;
            break;
        }
        default:
            break;
    }
    vector<unsigned> left_vertices_tmp(n_vertices[0]);
    vector<unsigned> right_vertices_tmp(n_vertices[1]);
    unsigned left_id = 0;
    unsigned right_id = 0;
    for(unsigned i = 0; i < all_vertices.size(); i++){
        if(all_vertices[i] < n_vertices[0]){
            left_vertices_tmp[left_id++] = all_vertices[i];
        }
        else{
            right_vertices_tmp[right_id++] = all_vertices[i];
        }
    }
    
    all_vertices.clear();
    if(anchor_left){
        all_vertices = left_vertices_tmp;
        all_vertices.insert(all_vertices.end(), right_vertices_tmp.begin(), right_vertices_tmp.end());
    }
    else{
        all_vertices = right_vertices_tmp;
        all_vertices.insert(all_vertices.end(), left_vertices_tmp.begin(), left_vertices_tmp.end());
    }
}

void SpecialBigraph::listing_cliques(){
    unsigned start_idx, max_d, max_two_hop_d, e, two_hop_e;
    unsigned *tmp_two_hop_d;
    
    //if(p > q){}
    //cout << "p=" << p << endl;
    //cout << "q=" << q << endl;
    n = n_vertices[0];
    
    d = new unsigned[n];
    tmp_two_hop_d = new unsigned[n];
    
    e = two_hop_e = 0;
    for(unsigned i = 0; i < n; i++){
        d[i] = adj_vec[i].size();
        e += d[i];
    }
    for(unsigned i = 0; i < n; i++){
        tmp_two_hop_d[i] = two_hop_adj_size[i];
        two_hop_e += tmp_two_hop_d[i];
    }
    
    
    cd = new unsigned[n+1];
    adj = new unsigned[e];
    two_hop_d = new unsigned*[p+1];
    two_hop_cd = new unsigned[n+1];
    two_hop_adj = new unsigned[two_hop_e];
    result_count = 0;
    
    
    cd[0] = 0;
    two_hop_cd[0] = 0;
    max_d = 0;
    max_two_hop_d = 0;
    
    for(unsigned i = 1; i < n+1; i++){
        cd[i] = cd[i-1] + d[i-1];
        max_d = (max_d > d[i-1])?max_d:d[i-1];
        two_hop_cd[i] = two_hop_cd[i-1] + tmp_two_hop_d[i-1];
        max_two_hop_d = (max_two_hop_d > tmp_two_hop_d[i-1])?max_two_hop_d:tmp_two_hop_d[i-1];
    }
        
    two_hop_core = max_two_hop_d;
    core = max_d;
    
    for(unsigned i = 0; i < n; i++){
        for(unsigned j = 0; j < adj_vec[i].size(); j++){
            adj[cd[i] + j] = adj_vec[i][j];
        }
    }
    
    //for(unsigned i = 0; i < e; i++){
    //    cout << adj[i] << " ";
    //}
    //cout << endl;
    
    for(unsigned i = 0; i < n; i++){
        for(unsigned j = 0; j < two_hop_adj_size[i]; j++){
            two_hop_adj[two_hop_cd[i] + j] = two_hop_adj_vec[i][j];
        }
    }
    
    //for(unsigned i = 0; i < two_hop_e; i++){
    //    cout << two_hop_adj[i] << " ";
    //}
    //cout << endl;
    
    #ifndef FREESUB
    printf("construct graph done\n"); ;
    #endif

    pqclique_main(p);
}

SpecialBigraph::Subgraph* SpecialBigraph::allocsub(unsigned l){
    
    Subgraph* sg = new Subgraph();

    #ifdef FREESUB
    //新加入的释放机制
    sg->freel = l;
    #endif


    sg->n = new unsigned[l+1]();
    sg->two_hop_d = new unsigned*[l+1];
    sg->two_hop_adj = new unsigned[two_hop_core * two_hop_core];
    sg->lab = new unsigned[two_hop_core];
    sg->nodes = new unsigned*[l+1];
    sg->op_vertices = new unsigned*[l+1];
    sg->op_size = new unsigned[l+1];
    sg->two_hop_core = two_hop_core;
    
    for(unsigned i = 1; i <= l; i++){
        sg->two_hop_d[i] = new unsigned[two_hop_core];
        sg->nodes[i] = new unsigned[two_hop_core];
        sg->op_vertices[i] = new unsigned[core];
    }
    
    //sg->new_id = new unsigned[n];
    sg->old_id = new unsigned[two_hop_core];
    
    return sg;
}

void SpecialBigraph::mksub(unsigned u, Subgraph* sg, unsigned l){
    unsigned i,j,k,v,w;
    
    #ifndef FREESUB
    static unsigned *new_id=NULL;
    #pragma omp threadprivate(new_id)
    if(new_id == NULL){
        new_id = new unsigned[n];
        for(unsigned i = 0; i < n; i++){
            new_id[i] = -1;
        }
    }
    #endif
    
    #ifdef FREESUB
    unsigned* new_id = new unsigned[n]; // 分配新的内存
    for(unsigned i = 0; i < n; i++) {
        new_id[i] = -1; // 初始化为 -1
    }
    #endif


    for(i = 0; i < sg->n[l-1]; i++){
        sg->lab[i] = 0;
    }
    
    for(i = 0; i <= l; i++){
        sg->op_size[i] = 0;
    }
        
    j = 0;
    for(i = two_hop_cd[u]; i < two_hop_cd[u+1]; i++){
        v = two_hop_adj[i];
        new_id[v] = j;
        // assert(j<1000);
        sg->old_id[j] = v;
        sg->lab[j] = l-1;
        sg->nodes[l-1][j] = j;
        sg->two_hop_d[l-1][j] = 0;
        j++;
    }
    
    // assert(j<1000);
    sg->n[l-1] = j;
        
    for(i = 0; i < sg->n[l-1]; i++){
        v = sg->old_id[i];
        for(k = two_hop_cd[v]; k < two_hop_cd[v+1]; k++){
            w = two_hop_adj[k];
            j = new_id[w];
            if(j != -1){
                // assert(j<1000) ;//
                sg->two_hop_adj[sg->two_hop_core * i + sg->two_hop_d[l-1][i]++] = j;
            }
        }
    }
    
    for(i = two_hop_cd[u]; i < two_hop_cd[u+1]; i++){
        v = two_hop_adj[i];
        new_id[v] = -1;
    }
    
    for(j = cd[u]; j < cd[u] + d[u]; j++){
        sg->op_vertices[l][sg->op_size[l]++] = adj[j];
    }

    #ifdef FREESUB
    delete [] new_id;
    #endif


}

void SpecialBigraph::pqclique_thread(unsigned l, Subgraph* sg, lint* rc){
    //cout << "sg->n[" << l << "]=" << sg->n[l] << endl;

    unsigned a,i,j,k,end,u,v,w;
    if(l == 2){
        #ifdef COUNT_ONLY
        for(i = 0; i < sg->n[2]; i++){
            u = sg->nodes[2][i];
            sg->op_size[2] = 0;
            end = cd[sg->old_id[u]] + d[sg->old_id[u]];
            
            j = 0;
            k = cd[sg->old_id[u]];
            while(j < sg->op_size[3] && k < cd[sg->old_id[u]] + d[sg->old_id[u]]){
                if(sg->op_vertices[3][j] == adj[k]){
                    sg->op_vertices[2][sg->op_size[2]++] = adj[k];
                    j++;
                    k++;
                }
                else if(sg->op_vertices[3][j] < adj[k]){
                    j++;
                }
                else{
                    k++;
                }
            }

            if(sg->op_size[l] < q){
                continue;
            }

            end = u * sg->two_hop_core + sg->two_hop_d[2][u];
            for(a = u * sg->two_hop_core; a < end; a++){
                v = sg->two_hop_adj[a];
                unsigned right_count = 0;
                j = 0;
                k = cd[sg->old_id[v]];
                while(j < sg->op_size[2] && k < cd[sg->old_id[v]] + d[sg->old_id[v]]){
                    if(sg->op_vertices[2][j] == adj[k]){
                        right_count++;
                        j++;
                        k++;
                    }
                    else if(sg->op_vertices[2][j] < adj[k]){
                        j++;
                    }
                    else{
                        k++;
                    }
                }
                *rc += Tools::choose(right_count, q);
            }
        }
        #else

        #endif
        //cout << "rc=" << *rc << endl;
        return;
    }
    else if(l == 1){
        #ifdef COUNT_ONLY
        for(i = 0; i < sg->n[1]; i++){
            unsigned right_count = 0;
            u = sg->nodes[1][i];
            j = 0;
            k = cd[sg->old_id[u]];
            while(j < sg->op_size[2] && k < cd[sg->old_id[u]] + d[sg->old_id[u]]){
                if(sg->op_vertices[2][j] == adj[k]){
                    right_count++;
                    j++;
                    k++;
                }
                else if(sg->op_vertices[2][j] < adj[k]){
                    j++;
                }
                else{
                    k++;
                }
            }
            *rc += Tools::choose(right_count, q);
        }
        #else

        #endif
        return;
    }
    for(i = 0; i < sg->n[l]; i++){
        u = sg->nodes[l][i];
        //cout << "u[" << l << "]=" << sg->old_id[u] << endl;
        
        // to handle opposite vertices;
        sg->op_size[l] = 0;
        j = 0;
        k = cd[sg->old_id[u]];
        while(j < sg->op_size[l+1] && k < cd[sg->old_id[u]] + d[sg->old_id[u]]){
            if(sg->op_vertices[l+1][j] == adj[k]){
                sg->op_vertices[l][sg->op_size[l]++] = adj[k];
                j++;
                k++;
            }
            else if(sg->op_vertices[l+1][j] < adj[k]){
                j++;
            }
            else{
                k++;
            }
        }
        
        if(sg->op_size[l] < q){
            continue;
        }

        sg->n[l-1] = 0;
        end = u * sg->two_hop_core + sg->two_hop_d[l][u];
        for(j = u * sg->two_hop_core; j < end; j++){
            v = sg->two_hop_adj[j];     //cout << "v=" << v << endl;
            if(sg->lab[v] == l){
                sg->lab[v] = l-1;
                sg->nodes[l-1][sg->n[l-1]++] = v;
                sg->two_hop_d[l-1][v] = 0;
            }
        }
        for(j = 0; j < sg->n[l-1]; j++){
            v = sg->nodes[l-1][j];
            end = sg->two_hop_core*v + sg->two_hop_d[l][v];
            for(k = sg->two_hop_core*v; k < end; k++){
                w = sg->two_hop_adj[k];
                if(sg->lab[w] == l-1){
                    sg->two_hop_d[l-1][v]++;
                }
                else{
                    sg->two_hop_adj[k--] = sg->two_hop_adj[--end];
                    sg->two_hop_adj[end] = w;
                }
            }
        }

        pqclique_thread(l-1, sg, rc);

        for(j = 0; j < sg->n[l-1]; j++){
            v = sg->nodes[l-1][j];
            sg->lab[v] = l;
        }
    }
}

void SpecialBigraph::pqclique_main(unsigned l){
    unsigned u;
    lint rc = 0;
    Subgraph *sg;
    #pragma omp parallel private(sg, u) reduction(+:rc)
    {
        sg = allocsub(l);
        #pragma omp for schedule(dynamic) nowait
        for(u = 0; u < n; u++){
            //cout << "root = " << u << endl;
            mksub(u, sg, l);
            pqclique_thread(l-1, sg, &rc);
            
            //cout << "result_count = " << rc << endl;
        }

        #ifdef FREESUB
        //新加入的释放机制
        delete sg; 
        #endif
        


    }
    result_count = rc;
}

SpecialBigraph::~SpecialBigraph(){


    //read_graph
    delete[] deg;
    vertices[0].clear();
    vertices[1].clear();
    vertices_in_left.clear();
    vertices_in_right.clear();
    for(unsigned i = 0; i < adj_vec.size(); i++){
        adj_vec[i].clear();
    }
    adj_vec.clear();
    edges.clear();
    list_of_edges.clear();

    if(if_prepare){
        //collect2
        for(unsigned i = 0; i < n_vertices[0]; i++){
            delete[] two_hop_adj_vec[i];
        }
        delete[] two_hop_adj_vec;
        delete[] two_hop_adj_size;
        delete[] two_hop_adj_maxsize;

        //mksub
        delete[] d;
        delete[] cd;
        delete[] adj;
        // for(unsigned i = 1; i <= p; i++){
        //     //delete[] two_hop_d[i];
        // }

        //list_clique
        delete[] two_hop_d;
        delete[] two_hop_cd;
        delete[] two_hop_adj;
    }


}


// the following are tool functions

void SpecialBigraph::print(vector<unsigned> vec){
    for(vector<unsigned>::iterator it = vec.begin(); it != vec.end(); it++){
        printf("%u ",*it);
    }
    printf("\n");
}

void SpecialBigraph::print_adj(){
    for(unsigned i = 0; i < num_vertices; i++){
        cout << "adj[" << i << "] : ";
        printf("adj[%u]:",i);
        for(unsigned j = 0; j < deg[i]; j++){
            printf("%d ",adj_vec[i][j]);
        }
        printf("\n");
    }
}

void SpecialBigraph::print_adj(unsigned i){
    for(unsigned j = 0; j < deg[i]; j++){
        printf("%d ",adj_vec[i][j]);
    }
    printf("\n");
}

void SpecialBigraph::print_two_hop_adj(){
    /*unsigned maxDeg = 0;
    
    for(unsigned i = 0; i < num_vertices; i++){
        maxDeg = maxDeg >= two_hop_adj_vec[i].size() ? maxDeg : two_hop_adj_vec[i].size();
        
        cout << "two_hop_adj[" << i << "] : ";
        for(unsigned j = 0; j < two_hop_adj_vec[i].size(); j++){
            cout << two_hop_adj_vec[i][j] << " ";
        }
        cout << endl;
    }
    cout << "max two hop degree : " << maxDeg << endl;*/
    
    unsigned maxDeg = 0;
    
    for(unsigned i = 0; i < n_vertices[0]; i++){
        maxDeg = maxDeg >= two_hop_adj_size[i] ? maxDeg : two_hop_adj_size[i];
        
        printf("two_hop_adj[%d] : ",i);
        for(unsigned j = 0; j < two_hop_adj_size[i]; j++){
            printf("%d ",two_hop_adj_vec[i][j]);
        }
        printf("\n");
    }
    printf("max two hop degree : %d\n",maxDeg);

}

void SpecialBigraph::print_edges(){
    cout << "#vertices : " << num_vertices <<endl;
    cout << "#edges : " << num_edges << endl;
    for(vector< pair<unsigned, unsigned> >::iterator it = list_of_edges.begin(); it != list_of_edges.end(); it++){
        cout << it->first << " " << it->second << endl;
    }
}

void SpecialBigraph::print_deg(){
    sort(deg, deg+num_vertices);
    for(unsigned i = num_vertices - 1; i >= num_vertices - 1000; i--){
        cout << deg[i] << endl;
    }
}

void SpecialBigraph::print_bclique(vector<unsigned> &left_vertices, vector<unsigned> &right_vertices){
    for(unsigned i = 0; i < left_vertices.size(); i++){
        printf("%d ",left_vertices[i]);
    }
    cout << "| ";
    for(unsigned i = 0; i < right_vertices.size(); i++){
        printf("%d ",right_vertices[i]);
    }
    printf("\n");
}

void SpecialBigraph::print_results(){
    #ifndef FREESUB
    cout << "Total # results : " << result_count << endl;
    #endif
    #ifndef COUNT_ONLY
    for(vector< pair<vector<unsigned>, vector<unsigned> > >::iterator it = bcliques.begin(); it != bcliques.end(); it++){
        print_bclique(it->first, it->second);
    }
    #endif
}



///// 接口

SpecialBigraph::SpecialBigraph(unsigned p_value, unsigned q_value, unsigned priority_strategy){

    clear_everything();
    p = p_value;
    q = q_value;
    priority = priority_strategy;

}//这里封装成一个vector



/*
这里改一个接口出来
*/

//测试接口
long long  getCountResult(string path , int p , int q , int strategy, int threadN){

    double beg_time, end_time, elapsed_time;
    
    
    SpecialBigraph *sbgraph = new SpecialBigraph(path,p,q,strategy);
    
    omp_set_num_threads(threadN);
    
    sbgraph->read_graph();

    
    // beg_time = time(NULL);
    if(!sbgraph->prepare_graph())
        return 0;
    // end_time = time(NULL);
    // elapsed_time = (end_time - beg_time);
    // cout << "Construct graph: using " << elapsed_time << " seconds." << endl;
    
    //return 0;
    
    // beg_time = time(NULL);

    long long re = 0;

    if(p==1&&q==1){
        re = sbgraph->num_edges; 
    }
    if(p==1){

    }else if (q==1) {
    
    }else {
        sbgraph->listing_cliques();
        re=sbgraph->result_count;
    }
   
    


    // end_time = time(NULL);
    // elapsed_time = (end_time - beg_time);
    // cout << "Computing biclique: using " << elapsed_time << " seconds." << endl;
    
    //cout << "Cost in finding adjs : " << time_counter1 << " seconds." << endl;
    // sbgraph->print_results();

    


    delete sbgraph;
    return re;
    

}

unsigned long long binomialCoefficient(int n, int k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;

    // 确保 k 是较小的那个
    if (k > n - k) k = n - k;

    unsigned long long result = 1;
    for (int i = 0; i < k; ++i) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

//程序接口
long long  getCountResult(vector<pair<uint, uint>> & edges  , int p , int q , int strategy, int threadN){

    double beg_time, end_time, elapsed_time;
    
    SpecialBigraph *sbgraph = new SpecialBigraph(p,q,strategy);
    
    omp_set_num_threads(threadN);
    omp_set_dynamic(0);
    
    sbgraph->read_graph(edges);//新的读取策略
    
    // beg_time = time(NULL);

    // end_time = time(NULL);
    // elapsed_time = (end_time - beg_time);
    // cout << "Construct graph: using " << elapsed_time << " seconds." << endl;
    
    //return 0;
    
    // beg_time = time(NULL);

    long long re = 0;
    if(p==0||q==0){
        re = 0 ; 
    }else if(p==1&&q==1){
        re = sbgraph->num_edges; 
    }
    if(p==1){

        //这里开始遍历        
        for(int i = 0 ; i < sbgraph->n_vertices[0] ; i++){
            if(sbgraph->deg[i] >= q){
                re+=binomialCoefficient(sbgraph->deg[i],q);
            }
        } 

    }else if (q==1) {
        
        //这里开始遍历        
        for(int i = sbgraph->n_vertices[0] ; i < sbgraph->num_vertices ; i++){
            if(sbgraph->deg[i] >= p){
                //总数取q加到re里面
                re+=binomialCoefficient(sbgraph->deg[i],p);
            }
        } 

    }else {
        if(!sbgraph->prepare_graph()) return 0;
        sbgraph->listing_cliques();
        re=sbgraph->result_count;
    }

    // end_time = time(NULL);
    // elapsed_time = (end_time - beg_time);
    // cout << "Computing biclique: using " << elapsed_time << " seconds." << endl;
    
    //cout << "Cost in finding adjs : " << time_counter1 << " seconds." << endl;


    // sbgraph->print_results();
        // 动态分配内存
        
    delete sbgraph;

    return re;
    
}

void SpecialBigraph::read_graph(vector<pair<uint, uint>> & edgeG){

    for(auto a : edgeG){
        add_edge(a.first, a.second);
    }

    vertices[0].clear();
    vertices[1].clear();
    
    n_vertices[0] = vertices_in_left.size();
    n_vertices[1] = vertices_in_right.size();
    num_vertices = n_vertices[0] + n_vertices[1];
    num_edges = edges.size();
    deg = new unsigned[num_vertices]();
    adj_vec.resize(num_vertices, vector<unsigned>());
    
    largest_index_in_partition[0] = 0;
    largest_index_in_partition[1] = n_vertices[0];
    for(set < pair <unsigned, unsigned> >::iterator it = edges.begin(); it != edges.end(); it++){
        unsigned A = it->first;
        unsigned B = it->second;
        get_index(A, 0);
        get_index(B, 1);
        adj_vec[A].emplace_back(B);
        adj_vec[B].emplace_back(A);
        list_of_edges.emplace_back(make_pair(A, B));
        deg[A]++;
        deg[B]++;
    }
    
    #ifdef IS_DEBUGGING
    //print_adj();
    
    //print_edges();
    #endif

    #ifndef FREESUB 
    printf("end of read graph\n");
    printf("#vertices = %u, #left_vertices = %u, #right_vertices = %u",num_vertices,n_vertices[0],n_vertices[1]);
    printf("#edges = %d\n",num_edges);    
    #endif

}


// int main(int argc, char *argv[]){
//     double beg_time, end_time, elapsed_time;
    
//     printf( "argc=%d\n", argc );
//     for( int i = 0; i < argc; ++i )
//         printf( "argv[%d]=%s\n", i, argv[i] );
    
//     SpecialBigraph *sbgraph = new SpecialBigraph(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
    
//     omp_set_num_threads(atoi(argv[5]));
    
//     beg_time = time(NULL);
//     sbgraph->read_graph();
//     end_time = time(NULL);
//     elapsed_time = (end_time - beg_time);
//     cout << "Read graph: using " << elapsed_time << " seconds." << endl;
    
//     beg_time = time(NULL);
//     if(!sbgraph->prepare_graph())
//         return 0;
//     end_time = time(NULL);
//     elapsed_time = (end_time - beg_time);
//     cout << "Construct graph: using " << elapsed_time << " seconds." << endl;
    
//     //return 0;
    
//     beg_time = time(NULL);
//     sbgraph->listing_cliques();
//     end_time = time(NULL);
//     elapsed_time = (end_time - beg_time);
//     cout << "Computing biclique: using " << elapsed_time << " seconds." << endl;
    
//     //cout << "Cost in finding adjs : " << time_counter1 << " seconds." << endl;
//     sbgraph->print_results();
//     delete sbgraph;
//     return 0;
// }




