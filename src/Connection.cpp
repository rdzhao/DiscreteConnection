#include "Connection.h"

Connection::Connection(){

}

Connection::~Connection(){

}

Polyhedron& Connection::Mesh(){
    return m_mesh;
}

int Connection::VertexIdx(Vertex_iterator vi){
    return m_vertex_idx[vi];
}

int Connection::EdgeIdx(Edge_iterator ei){
    return m_edge_idx[ei];
}

int Connection::HalfedgeIdx(Halfedge_iterator hei){
    return m_halfedge_idx[hei];
}

int Connection::FacetIdx(Facet_iterator fi){
    return m_facet_idx[fi];
}

double Connection::RhoE(Halfedge_iterator hei){
    return m_rho_e[hei];
}

double Connection::RhoVE(Halfedge_iterator hei){
    return m_rho_ve[hei];
}

double Connection::RhoVT(Halfedge_iterator hei){
    return m_rho_vt[hei];
}

void Connection::read(std::string mesh_fn){
    readMesh(mesh_fn);    
}

void Connection::preprocessing(){
    indexElements();
    computeAverageEdgeLength();

    //testBorderEdge();
    //testInnerProduct();
    //testCotan();
    //testAngleDefect();
    //testOrientedAngle();

    //totalAngleDefect();
}

void Connection::compute(){
    std::cout<<"Assembling Linear System ..."<<std::endl;
    assembleLinearSystem();
    std::cout<<"Solving Linear System ..."<<std::endl;
    solveLinearSystem();
    std::cout<<"Computing Field By Connection..."<<std::endl;
    computeFieldByConnection();
    saveVertexSpanningTree();
    std::cout<<"Converting to Vertex-based Vector Field ..."<<std::endl;
    computeVertexBasedVF();

    //checkDT();
    //checkDE();
    //checkLCCInducedDefect();
    //checkDefect();
    
    //testHalfedgeRotation();
    //testFacetRotation();
}

void Connection::write(){
    std::cout<<"Writing ... "<<std::endl;
    //writeMeshToOFF();
    writeMeshToVTK();
    writeMeshToVTKWithField();

    //saveFrames();
    //saveVertexToEdgeFrames();
    //saveVertexToTriangleFrames();
    //std::cout<<"LCF ..."<<std::endl;
    //computeOptimalHalfedgeBasedConnection();
    //testLeviCivitaConnection();
    //saveFacetSpanningTree();

    //saveVertexOneRingRotation();
    //saveVertexOneRingUniversalTest();
}

void Connection::readMesh(std::string mesh_fn){
    std::ifstream file(mesh_fn.c_str());
    file >> m_mesh;
}

void Connection::indexElements(){
    int numv = 0;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        m_vertex_idx[vi] = numv++;
    }

    int nume = 0;
    for(Edge_iterator ei = m_mesh.edges_begin(); ei !=m_mesh.edges_end(); ++ei){
        m_edge_idx[ei] = nume++;
    }

    int numhe = 0;
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        m_halfedge_idx[hei] = numhe++;
    }

    int numf = 0;
    for(Facet_iterator fi =m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        m_facet_idx[fi] = numf++;
    }

    std::cout<<"Num of Vertex: "<<m_mesh.size_of_vertices()<<std::endl;
    std::cout<<"Num of Halfedge: "<<m_mesh.size_of_halfedges()<<std::endl;
    std::cout<<"Num of Facet: "<<m_mesh.size_of_facets()<<std::endl;
}

void Connection::computeAverageEdgeLength(){
    m_ael = 0;
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        m_ael += length(ei);
    }
    m_ael /= m_mesh.size_of_halfedges()/2;
}

void Connection::assembleLinearSystem(){
    std::cout<<"Computing Levi-Civita Connection ..."<<std::endl;
    computeLeviCivitaConnection();
    std::cout<<"Determining Integer Jump ..."<<std::endl;
    determineIntegerJump();
    
    std::cout<<"Assembling Matrix ..."<<std::endl;
    assembleEnergyMatrix();
    std::cout<<"Assembling Vector ..."<<std::endl;
    assembleEnergyVector();
}

void Connection::solveLinearSystem(){
    
    fixGauge();
    saveEnergyMatrix();

    Eigen::SimplicialLDLT<SparseMatrix> solver;
    solver.compute(m_A);
    if(solver.info()!=Eigen::Success) {
        std::cout<<"Compute Fail ..."<<std::endl;
        std::abort();
    }
    VectorXd x = solver.solve(m_b/2);
    if(solver.info()!=Eigen::Success) {
        std::cout<<"Solve Fail ..."<<std::endl;
        std::abort();
    }

    //std::cout<<x<<std::endl;
    saveVector(x);

    // store the results
    int block_1 = 3*m_mesh.size_of_facets();
    int block_12 = block_1 + m_mesh.size_of_halfedges();
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        m_rho_vt[hei] = x(m_halfedge_idx[hei]);
        m_rho_ve[hei] = x(block_1+m_halfedge_idx[hei]);
    }
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        m_rho_e[ei] = x(block_12+m_edge_idx[ei]);
    }
}

void Connection::computeFieldByConnection(){
    // init spanning tree
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        m_inst[ei] = false;
    }

    std::map<Vertex_iterator, bool> explored;
    std::map<Vertex_iterator, bool> visited;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        explored[vi] = false;
        visited[vi] = false;
    }

    std::list<Vertex_iterator> toBeExplored;
    Vertex_iterator bvi = m_mesh.vertices_begin();
    //for(int k=0; k<30; ++k){
    //    ++bvi;
    //}
    //std::cout<<m_vertex_idx[bvi]<<std::endl;
    toBeExplored.push_back(bvi);
    m_field[bvi] = 0;
    visited[bvi] = true;

    //int iter = 0;
    while(!toBeExplored.empty()){
        Vertex_iterator vi = toBeExplored.front();
        toBeExplored.pop_front();
        explored[vi] = true;

        Halfedge_around_vertex_circulator hecir = vi->vertex_begin();
        do{
            Halfedge_iterator thei = hecir;
            if(!visited[thei->opposite()->vertex()]){
                toBeExplored.push_back(thei->opposite()->vertex());

                if(!visited[thei->opposite()->vertex()]){
                    if(isEdge(thei)){
                        m_field[thei->opposite()->vertex()] = m_field[vi] - (- m_rho_e[thei] + m_q[thei->opposite()]*2*PI);
                        m_inst[thei] = true;
                    }
                    else{
                        m_field[thei->opposite()->vertex()] = m_field[vi] - m_rho_e[thei->opposite()];
                        m_inst[thei->opposite()] = true;
                    }
                }
                visited[thei->opposite()->vertex()] = true;
            }
            ++hecir;
            //++iter;
        }while(hecir != vi->vertex_begin());
    }

    //for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
    //    if(counter[vi] != 1){
    //        std::cout<<"WTF is this !!! "<<m_vertex_idx[vi]<<std::endl;
    //    }
    //}
}

void Connection::computeVertexBasedVF(){
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        Vector_3 f1, f2;
        frame(f1, f2, vi);

        Vector_3 v = std::cos(m_field[vi])*f1 + std::sin(m_field[vi])*f2;
        m_vf[vi] = v;
    }
}

void Connection::writeMeshToOFF(){
    std::ofstream out("mesh.off");

    // info
    out<<"COFF"<<std::endl;
    out<<m_mesh.size_of_vertices()<<" "<<m_mesh.size_of_facets()<<" 0"<<std::endl;

    // vertices
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi!=m_mesh.vertices_end(); ++vi){
        out<<vi->point().x()<<" "
            <<vi->point().y()<<" "
            <<vi->point().z()<<std::endl;
    }

    // facets
    for(Facet_iterator fi = m_mesh.facets_begin(); fi!=m_mesh.facets_end(); ++fi){
        out<<"3 "
            <<m_vertex_idx[fi->halfedge()->vertex()]<<" "
            <<m_vertex_idx[fi->halfedge()->next()->vertex()]<<" "
            <<m_vertex_idx[fi->halfedge()->next()->next()->vertex()]<<std::endl;
    }

    out.close();
}

void Connection::writeMeshToVTK(){
    std::ofstream out("mesh.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Mesh with vertex based vector field"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    out<<"POINTS "<<m_mesh.size_of_vertices()<<" double"<<std::endl;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        out<<vi->point().x()<<" "
            <<vi->point().y()<<" "
            <<vi->point().z()<<std::endl;
    }

    out<<"CELLS "<<m_mesh.size_of_facets()<<" "<<4*m_mesh.size_of_facets()<<std::endl;
    for(Facet_iterator fi = m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        out<<"3 "
            <<m_vertex_idx[fi->halfedge()->vertex()]<<" "
            <<m_vertex_idx[fi->halfedge()->next()->vertex()]<<" "
            <<m_vertex_idx[fi->halfedge()->next()->next()->vertex()]<<std::endl;
    }

    out<<"CELL_TYPES "<<m_mesh.size_of_facets()<<std::endl;
    for(Facet_iterator fi = m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        out<<"5"<<std::endl;
    }

    out.close();
}

void Connection::writeMeshToVTKWithField(){
    std::ofstream out("mesh_vf.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Mesh with vertex based vector field"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    out<<"POINTS "<<m_mesh.size_of_vertices()<<" double"<<std::endl;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        out<<vi->point().x()<<" "
            <<vi->point().y()<<" "
            <<vi->point().z()<<std::endl;
    }

    out<<"CELLS "<<m_mesh.size_of_facets()<<" "<<4*m_mesh.size_of_facets()<<std::endl;
    for(Facet_iterator fi = m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        out<<"3 "
            <<m_vertex_idx[fi->halfedge()->vertex()]<<" "
            <<m_vertex_idx[fi->halfedge()->next()->vertex()]<<" "
            <<m_vertex_idx[fi->halfedge()->next()->next()->vertex()]<<std::endl;
    }

    out<<"CELL_TYPES "<<m_mesh.size_of_facets()<<std::endl;
    for(Facet_iterator fi = m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        out<<"5"<<std::endl;
    }

    out<<"POINT_DATA "<<m_mesh.size_of_vertices()<<std::endl;
    out<<"VECTORS U double"<<std::endl;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        out<<m_vf[vi].x()<<" "<<m_vf[vi].y()<<" "<<m_vf[vi].z()<<std::endl;
    }

    out.close();
}

void Connection::computeLeviCivitaConnection()
{
    // Levi-Civita connection is defined by oriented angle betwen [0, 2*PI).
    // If transfer from e to t, one should subtract the computed angle here as e->t.
    for(Halfedge_iterator hei=m_mesh.halfedges_begin(); hei!=m_mesh.halfedges_end(); ++hei){
        Halfedge_iterator fhei = hei->facet()->halfedge();
        m_levi_civita[hei] = orientedAngle(hevector(hei), hevector(fhei), normal(hei->facet()));
    }
}

void Connection::determineIntegerJump(){
    determine_n();
    determine_q();
}

void Connection::determine_n(){
    int l,m;
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        // determining integer jump l associated with
        // tail of hei -> opposite of hei
        if(hei->opposite() == hei->opposite()->facet()->halfedge()){
            l = 0;
        }
        else{
            l = 1;
        }

        // Take care: CGAL vertex->halfedge() is the halfedge that point to this vertex
        if(hei == hei->opposite()->vertex()->halfedge()->opposite()){
            m = 1;
        }
        else{
            m = 0;
        }

        m_n[hei] = m-l;
    }
}

void Connection::determine_q(){
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        if(m_halfedge_idx[hei]%2 == 0){ // is edge
            m_q[hei] = 0;
        }
        else{
            m_q[hei] = -1 - m_n[hei] - m_n[hei->opposite()];
        }
    }
}

void Connection::assembleEnergyMatrix(){
    std::vector<Triplet> triplets_Et;
    assembleRhovtRhovt_Et(triplets_Et);
    assembleRhovRhov_Et(triplets_Et);
    assembleRhovtRhov_Et(triplets_Et);
    
    m_SMEt.resize(2.5*m_mesh.size_of_halfedges(), 2.5*m_mesh.size_of_halfedges());
    m_SMEt.setFromTriplets(triplets_Et.begin(), triplets_Et.end());
    //isSymmetric(m_SMEt);
    //countNonzeroEntries(m_SMEt);

    // wrong !!! no cotan formula
    std::vector<Triplet> triplets_Ee;
    assembleRhovtRhovt_Ee(triplets_Ee);
    assembleRhoveRhove_Ee(triplets_Ee);
    assembleRhovtRhove_Ee(triplets_Ee);
    
    m_SMEe.resize(2.5*m_mesh.size_of_halfedges(), 2.5*m_mesh.size_of_halfedges());
    m_SMEe.setFromTriplets(triplets_Ee.begin(), triplets_Ee.end());
    //isSymmetric(m_SMEe);

    m_A = m_SMEt + m_SMEe;

    // debug
    //saveEnergyMatrix();
}

void Connection::assembleEnergyVector(){
    m_b.setZero(2.5*m_mesh.size_of_halfedges());
    assembleVectorFromEt(m_b);
    assembleVectorFromEe(m_b);

    //std::cout<<m_b<<std::endl;
}

void Connection::assembleRhovtRhovt_Et(std::vector<Triplet>& triplets){
    // halfedge storage order
    // each halfedge represents the tail vertex to associated facet
    // diagonal term
    for(Facet_iterator fi = m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        Halfedge_iterator hei = fi->halfedge();
        do{
            triplets.push_back(Triplet(m_halfedge_idx[hei], m_halfedge_idx[hei], 
                ipSame1F(hei)
                + ipSame1F(hei->prev())
                - 2*ipDiff1F(hei->prev())
                ));   
            hei = hei->next();
        }while(hei!=fi->halfedge());
    }

    // off diagonal term
    for(Facet_iterator fi = m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        Halfedge_iterator hei = fi->halfedge();
        do{
            triplets.push_back(Triplet(m_halfedge_idx[hei], m_halfedge_idx[hei->next()],
                            - ipSame1F(hei)
                            + ipDiff1F(hei)
                            + ipDiff1F(hei->prev())
                            - ipDiff1F(hei->next())
                            ));
            triplets.push_back(Triplet(m_halfedge_idx[hei], m_halfedge_idx[hei->prev()],
                            - ipDiff1F(hei)
                            + ipDiff1F(hei->prev())
                            + ipDiff1F(hei->next())
                            - ipSame1F(hei->prev())
                            ));

            hei = hei->next();
        }while(hei!=fi->halfedge());
    }
}

void Connection::assembleRhovRhov_Et(std::vector<Triplet>& triplets){
    // edge storage order
    // each edge with its orientation as halfedge
    // diagonal
    // each self-multiplication will only appear in two incident triangles
    int block_1 = 3*m_mesh.size_of_facets();
    int block_12 = block_1 + m_mesh.size_of_halfedges();
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        triplets.push_back(Triplet(block_12+m_edge_idx[ei], block_12+m_edge_idx[ei], 
            ipSame1F(ei)
            + ipSame1F(ei->opposite())
            ));
    }

    // off diagonal
    // each edge will have interacton with 4 edges.
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        Edge_iterator nei = ei->next();
        int ns = 1; // because ei has secured positive sign
        if(m_halfedge_idx[nei]%2 == 1){
            nei = nei->opposite();
            ns = -1;
        }
        triplets.push_back(Triplet(block_12 + m_edge_idx[ei], block_12 + m_edge_idx[nei], ns*ipDiff1F(ei)));

        Edge_iterator pei = ei->prev();
        int ps = 1;
        if(m_halfedge_idx[pei]%2 == 1){
            pei = pei->opposite();
            ps = -1;
        }
        triplets.push_back(Triplet(block_12 + m_edge_idx[ei], block_12 + m_edge_idx[pei], ps*ipDiff1F(ei->prev())));

        Edge_iterator onei = ei->opposite()->next();
        int ons = -1; // because ei->opposite() has secured positive sign
        if(m_halfedge_idx[onei]%2 == 1){
            onei = onei->opposite();
            ons = 1;
        }
        triplets.push_back(Triplet(block_12 + m_edge_idx[ei], block_12 + m_edge_idx[onei], ons*ipDiff1F(ei->opposite())));

        Edge_iterator opei = ei->opposite()->prev();
        int ops = -1;
        if(m_halfedge_idx[opei]%2 == 1){
            opei = opei->opposite();
            ops = 1;
        }
        triplets.push_back(Triplet(block_12 + m_edge_idx[ei], block_12 + m_edge_idx[opei], ops*ipDiff1F(ei->opposite()->prev())));
    }
}

void Connection::assembleRhovtRhov_Et(std::vector<Triplet>& triplets){
    // iterator through all the halfedges as v->t
    // each v->t will have interactions will all the 3 edges in the corresponding triangle
    int block_1 = 3*m_mesh.size_of_facets();
    int block_12 = block_1 + m_mesh.size_of_halfedges();
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        // hei here represents a instance of v->t
        // need to first just halfedge whether an edge
        // first edge
        Halfedge_iterator thei1 = hei; // thei represents an instance of v->v
        int s1 = 1; // solely for determining whether v->v has negative sign, v->t sign is handled in triplets assembly.
        if(m_halfedge_idx[thei1]%2 != 0){
            thei1 = thei1->opposite();
            s1 = -1;
        }
        triplets.push_back(Triplet(m_halfedge_idx[hei], block_12 + m_edge_idx[thei1],
            -s1*ipSame1F(hei)
            +s1*ipDiff1F(hei->prev())
            ));
        // symmetric
        triplets.push_back(Triplet(block_12 + m_edge_idx[thei1], m_halfedge_idx[hei], 
            -s1*ipSame1F(hei)
            +s1*ipDiff1F(hei->prev())
            ));

        // second edge
        Halfedge_iterator thei2 = hei->next();
        int s2 = 1;
        if(m_halfedge_idx[thei2]%2 != 0){
            thei2 = thei2->opposite();
            s2 = -1;
        }
        triplets.push_back(Triplet(m_halfedge_idx[hei], block_12 + m_edge_idx[thei2],
            -s2*ipDiff1F(hei)
            +s2*ipDiff1F(hei->next())
            ));
        // symmetric
        triplets.push_back(Triplet(block_12 + m_edge_idx[thei2], m_halfedge_idx[hei],
            -s2*ipDiff1F(hei)
            +s2*ipDiff1F(hei->next())
            ));

        // thrid edge
        Halfedge_iterator thei3 = hei->prev();
        int s3 = 1;
        if(m_halfedge_idx[thei3]%2 != 0){
            thei3 = thei3->opposite();
            s3 = -1;
        }
        triplets.push_back(Triplet(m_halfedge_idx[hei], block_12 + m_edge_idx[thei3],
            -s3*ipDiff1F(hei->prev())
            +s3*ipSame1F(hei->prev())
        ));
        // symmetric
        triplets.push_back(Triplet(block_12 + m_edge_idx[thei3], m_halfedge_idx[hei],
            -s3*ipDiff1F(hei->prev())
            +s3*ipSame1F(hei->prev())
        ));
    }
}

void Connection::assembleRhovtRhovt_Ee(std::vector<Triplet>& triplets){
    // iterate through all the v->t represented as halfedge
    // v->t will only appear in the same triangle
    // thus only 2 edges that connect to this v in the same triangle will have contribution
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        // diagonal
        triplets.push_back(Triplet(m_halfedge_idx[hei], m_halfedge_idx[hei],
            ipSame0F(hei)*cotan(hei) // appear in edge to triangle using hei as representative (hei to hei->facet())
            + ipSame0F(hei->prev())*cotan(hei->prev())
            ));
    }

    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        // off diagonal
        triplets.push_back(Triplet(m_halfedge_idx[hei], m_halfedge_idx[hei->next()],
            ipDiff0F(hei)*cotan(hei)
            ));
        triplets.push_back(Triplet(m_halfedge_idx[hei], m_halfedge_idx[hei->prev()],
            ipDiff0F(hei->prev())*cotan(hei->prev())
            ));
    }
}

void Connection::assembleRhoveRhove_Ee(std::vector<Triplet>& triplets){
    // halfedge storage order
    // each halfedge represents tail vertex to associated halfedge
    // each halfedge has contributions from two adjacent triangles
    // diagonal term 
    int block_1 = 3*m_mesh.size_of_facets();
    int block_12 = block_1 + m_mesh.size_of_halfedges();
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        triplets.push_back(Triplet(block_1 + m_halfedge_idx[hei], block_1 + m_halfedge_idx[hei],
            ipSame0F(hei)*cotan(hei)
            + ipSame0F(hei->opposite())*cotan(hei->opposite())
            ));
    }

    // off diagonal term
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        triplets.push_back(Triplet(block_1 + m_halfedge_idx[hei], block_1 + m_halfedge_idx[hei->opposite()],
            ipDiff0F(hei)*cotan(hei)
            + ipDiff0F(hei->opposite())*cotan(hei->opposite())
            ));
    }
}

void Connection::assembleRhovtRhove_Ee(std::vector<Triplet>& triplets){
    // iterate through all the halfedges as v->t
    // each v->t will have interaction with two edges incident to this vertex
    int block_1 = 3*m_mesh.size_of_facets();
    int block_12 = block_1 + m_mesh.size_of_halfedges();
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        triplets.push_back(Triplet(m_halfedge_idx[hei], block_1 + m_halfedge_idx[hei],
            -ipSame0F(hei)*cotan(hei)
            ));
        triplets.push_back(Triplet(m_halfedge_idx[hei], block_1 + m_halfedge_idx[hei->opposite()],
            -ipDiff0F(hei)*cotan(hei)
            ));
        triplets.push_back(Triplet(m_halfedge_idx[hei], block_1 + m_halfedge_idx[hei->prev()],
            -ipDiff0F(hei->prev())*cotan(hei->prev())
            ));
        triplets.push_back(Triplet(m_halfedge_idx[hei], block_1 + m_halfedge_idx[hei->prev()->opposite()],
            -ipSame0F(hei->prev())*cotan(hei->prev())
            ));

        // symmetric 
        triplets.push_back(Triplet(block_1 + m_halfedge_idx[hei], m_halfedge_idx[hei], 
            -ipSame0F(hei)*cotan(hei)
            ));
        triplets.push_back(Triplet(block_1 + m_halfedge_idx[hei->opposite()], m_halfedge_idx[hei], 
            -ipDiff0F(hei)*cotan(hei)
            ));
        triplets.push_back(Triplet(block_1 + m_halfedge_idx[hei->prev()], m_halfedge_idx[hei], 
            -ipDiff0F(hei->prev())*cotan(hei->prev())
            ));
        triplets.push_back(Triplet(block_1 + m_halfedge_idx[hei->prev()->opposite()], m_halfedge_idx[hei], 
            -ipSame0F(hei->prev())*cotan(hei->prev())
            ));
    }
}

void Connection::assembleVectorFromEt(VectorXd& v){
    int block_1 = 3*m_mesh.size_of_facets();
    int block_12 = block_1 + m_mesh.size_of_halfedges();
    for(Facet_iterator fi = m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        Halfedge_iterator hei = fi->halfedge();
        do{
            // in front of v->t
            v(m_halfedge_idx[hei]) -= -2*m_q[hei]*2*PI*ipSame1F(hei);
            v(m_halfedge_idx[hei->next()]) -= 2*m_q[hei]*2*PI*ipSame1F(hei);

            v(m_halfedge_idx[hei->next()]) -= -m_q[hei]*2*PI*ipDiff1F(hei);
            v(m_halfedge_idx[hei->prev()]) -= m_q[hei]*2*PI*ipDiff1F(hei);
            v(m_halfedge_idx[hei]) -= -m_q[hei->next()]*2*PI*ipDiff1F(hei);
            v(m_halfedge_idx[hei->next()]) -= m_q[hei->next()]*2*PI*ipDiff1F(hei);

            v(m_halfedge_idx[hei->prev()]) -= -m_q[hei]*2*PI*ipDiff1F(hei->prev());
            v(m_halfedge_idx[hei]) -= m_q[hei]*2*PI*ipDiff1F(hei->prev());
            v(m_halfedge_idx[hei]) -= -m_q[hei->prev()]*2*PI*ipDiff1F(hei->prev());
            v(m_halfedge_idx[hei->next()]) -= m_q[hei->prev()]*2*PI*ipDiff1F(hei->prev());

            // in front of v->v
            //std::cout<<"Check: "<<isEdge(hei)<<" "<<m_q[hei]<<std::endl;
            int s_ij = 1;
            int s_jk = 1;
            int s_ki = 1;
            if(m_halfedge_idx[hei]%2 == 1){
                s_ij = -1;
            }
            if(m_halfedge_idx[hei->next()]%2 == 1){
                s_jk = -1;
            }
            if(m_halfedge_idx[hei->prev()]%2 == 1){
                s_ki = -1;
            }
            
            v(block_12 + edgeIdxFromHalfedge(hei)) -= s_ij*2*m_q[hei]*2*PI*ipSame1F(hei);

            v(block_12 + edgeIdxFromHalfedge(hei->next())) -= s_jk*m_q[hei]*2*PI*ipDiff1F(hei);
            v(block_12 + edgeIdxFromHalfedge(hei)) -= s_ij*m_q[hei->next()]*2*PI*ipDiff1F(hei);

            v(block_12 + edgeIdxFromHalfedge(hei->prev())) -= s_ki*m_q[hei]*2*PI*ipDiff1F(hei->prev());
            v(block_12 + edgeIdxFromHalfedge(hei)) -= s_ij*m_q[hei->prev()]*2*PI*ipDiff1F(hei->prev());

            hei=hei->next();
        }while(hei!=fi->halfedge());
    }
}

void Connection::assembleVectorFromEe(VectorXd& v){
    // iterator through all the halfedges as e_ij->t_ijk
    int block_1 = 3*m_mesh.size_of_facets();
    int block_12 = block_1 + m_mesh.size_of_halfedges();
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        // from quadratic term
        v(m_halfedge_idx[hei]) -= -2*(0.5+m_n[hei->opposite()])*2*PI*ipDiff0F(hei)*cotan(hei);
        v(block_1 + m_halfedge_idx[hei]) -= 2*(0.5+m_n[hei->opposite()])*2*PI*ipDiff0F(hei)*cotan(hei);
        
        v(m_halfedge_idx[hei->next()]) -= -2*(0.5+m_n[hei->opposite()])*2*PI*ipSame0F(hei)*cotan(hei);
        v(block_1 + m_halfedge_idx[hei->opposite()]) -= 2*(0.5+m_n[hei->opposite()])*2*PI*ipSame0F(hei)*cotan(hei);

        // from linear term
        v(m_halfedge_idx[hei]) -= -2*(0.5*length(hei))*m_levi_civita[hei]*cotan(hei);
        v(block_1+m_halfedge_idx[hei]) -= 2*(0.5*length(hei))*m_levi_civita[hei]*cotan(hei);

        v(m_halfedge_idx[hei->next()]) -= -2*(0.5*length(hei))*m_levi_civita[hei]*cotan(hei);
        v(block_1 + m_halfedge_idx[hei->opposite()]) -= 2*(0.5*length(hei))*m_levi_civita[hei]*cotan(hei);
    }

}

void Connection::fixGauge(){
    int block_1 = 3*m_mesh.size_of_facets();
    std::set<int> indices;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        indices.insert(block_1 + m_halfedge_idx[vi->halfedge()->opposite()]);
    }    

    for(int i=0; i<m_A.outerSize(); ++i){
        for(SparseMatrix::InnerIterator iter(m_A, i); iter; ++iter){
            if(indices.find(iter.row()) != indices.end() || indices.find(iter.col()) != indices.end()){
                if(iter.row() == iter.col()){
                    iter.valueRef() = 1;
                }   
                else{
                    iter.valueRef() = 0;
                }
            }
        }
    }

    for(std::set<int>::iterator iter = indices.begin(); iter!=indices.end(); ++iter){
        m_b(*iter) = 0;
    }
}

double Connection::ipSame0F(Halfedge_iterator hei){
    return length(hei)/3;
}

double Connection::ipDiff0F(Halfedge_iterator hei){
    return length(hei)/6;
}

double Connection::ipSame1F(Halfedge_iterator hei){
    double a = area(hei->facet());
    double val = (a/6)/(height(hei->next())*height(hei->next()))
                + (a/6)/(height(hei->prev())*height(hei->prev()))
                - 2*(a/12)/(height(hei->next())*height(hei->prev()));
    //std::cout<<"Look inside 1: "<<height(hei->next())*height(hei->next())<<std::endl;
    //std::cout<<"Look inside 2: "<<height(hei->prev())*height(hei->prev())<<std::endl;
    //std::cout<<"Look inside 3: "<<height(hei->next())*height(hei->prev())<<std::endl;
    //std::cout<<height(hei->next())<<" :: "<<height(hei->prev())<<std::endl; 
    return val;
}

double Connection::ipDiff1F(Halfedge_iterator hei){
    double a =area(hei->facet());
    double val = (a/12)/(height(hei->prev())*height(hei))
                - (a/6)/(height(hei->next())*height(hei))
                - (a/12)/(height(hei->prev())*height(hei->prev()))
                + (a/12)/(height(hei->next())*height(hei->prev()));
    return val;
}

double Connection::norm(Vector_3 v){
    return sqrt(v.squared_length());
}

Vector_3 Connection::cross(Vector_3 v1, Vector_3 v2){
    return CGAL::cross_product(v1,v2);
}

Point_3 Connection::baryCenter(Edge_iterator ei){
    return ei->opposite()->vertex()->point() + (ei->vertex()->point() - ei->opposite()->vertex()->point())/2;
}

Vector_3 Connection::normal(Edge_iterator ei){
    Vector_3 n = normal(ei->facet())*area(ei->facet()) + normal(ei->opposite()->facet())*area(ei->opposite()->facet());
    double ta  = area(ei->facet()) + area(ei->opposite()->facet());
    
    n /= ta;
    n /= norm(n);

    return n;
}

Vector_3 Connection::height(Halfedge_iterator hei){
    double a = area(hei->facet());
    double hm = 2*a/length(hei);
    
    Vector_3 n = normal(hei->facet());
    Vector_3 v = hei->vertex()->point() - hei->opposite()->vertex()->point();
    v/=sqrt(v.squared_length());
    Vector_3 hd = rotate(v, n, PI/2);

    return hm*hd;
}

double Connection::internalAngle(Halfedge_iterator hei){
    Vector_3 v1 = hei->opposite()->vertex()->point() - hei->vertex()->point();
    Vector_3 v2 = hei->next()->vertex()->point() - hei->vertex()->point();
    //std::cout<<std::acos(v1*v2/(norm(v1)*norm(v2)))<<std::endl;
    
    return std::acos(v1*v2/(norm(v1)*norm(v2)));
}

double Connection::length(Halfedge_iterator hei){
    return sqrt((hei->vertex()->point()-hei->opposite()->vertex()->point()).squared_length());
}

double Connection::area(Facet_iterator fi){
    Halfedge_iterator hei = fi->halfedge();
    Vector_3 v1 = hei->next()->vertex()->point()-hei->vertex()->point();
    Vector_3 v2 = hei->prev()->vertex()->point()-hei->vertex()->point();
    Vector_3 va = CGAL::cross_product(v1, v2);
    return sqrt(va.squared_length())/2;
}

Point_3 Connection::baryCenter(Facet_iterator fi){
    Vector_3 v1 = fi->halfedge()->vertex()->point() - CGAL::ORIGIN;
    Vector_3 v2 = fi->halfedge()->next()->vertex()->point() - CGAL::ORIGIN;
    Vector_3 v3 = fi->halfedge()->next()->next()->vertex()->point() - CGAL::ORIGIN;
    return CGAL::ORIGIN + (v1 + v2 + v3)/3;
}

Vector_3 Connection::normal(Facet_iterator fi){
    Halfedge_iterator hei = fi->halfedge();
    Vector_3 v1 = hei->next()->vertex()->point()-hei->vertex()->point();
    Vector_3 v2 = hei->prev()->vertex()->point()-hei->vertex()->point();
    Vector_3 va = CGAL::cross_product(v1, v2);
    return va / norm(va);
}

Vector_3 Connection::normal(Vertex_iterator vi){
    Halfedge_around_vertex_circulator hecir = vi->vertex_begin();
    double w = 0;
    Vector_3 n(0,0,0);
    do{
        n+=area(hecir->facet())*normal(hecir->facet());
        w+=area(hecir->facet());

        ++hecir;
    }while(hecir != vi->vertex_begin());

    n /= w;
    n /= norm(n);

    return n;
}

double Connection::defect(Vertex_iterator vi){
    Halfedge_iterator hei = vi->halfedge();
    double ta = 0;
    do{
        Vector_3 v1 = hei->opposite()->vertex()->point() - vi->point();
        Vector_3 v2 = hei->next()->vertex()->point() - vi->point();
        ta += angle(v1,v2);

        hei = hei->next()->opposite();
    }while(hei != vi->halfedge());

    return 2*PI - ta;
}

int Connection::valence(Vertex_iterator vi){
    Halfedge_iterator hei = vi->halfedge();
    int vl = 0;
    do{
        ++vl;

        hei = hei->next()->opposite();
    }while(hei != vi->halfedge());

    return vl;
}

double Connection::angle(Vector_3 v1, Vector_3 v2){
    double val = v1*v2/(norm(v1)*norm(v2));
    if(val > 1)
        val = 1;
    if(val < -1)
        val = -1;
    
    return std::acos(val);
}

double Connection::orientedAngle(Vector_3 v1, Vector_3 v2, Vector_3 n){        
    double y = cross(v1,v2)*n;
    double x = v1*v2;

    if(y>=0){
        return std::atan2(y,x);
    }
    else{
        return 2*PI+std::atan2(y,x);
    }
}

Point_3 Connection::pointOnHalfedge(Halfedge_iterator hei, double t){
    Vector_3 v1 = hei->vertex()->point() - CGAL::ORIGIN;
    Vector_3 v2 = hei->opposite()->vertex()->point() - CGAL::ORIGIN;
    return CGAL::ORIGIN + t*v1 + (1-t)*v2;
}

Point_3 Connection::pointOnFacet(Halfedge_iterator hei, double t){
    Vector_3 v1 = hei->vertex()->point() - CGAL::ORIGIN;
    Vector_3 v2 = hei->next()->vertex()->point() - CGAL::ORIGIN;
    Vector_3 v3 = hei->next()->next()->vertex()->point() - CGAL::ORIGIN;
    return CGAL::ORIGIN + t*v3 + (1-t)*(v1+v2)/2;
}

double Connection::angleVertexToEdge(Halfedge_iterator hei){
    // head vertex
    if(m_halfedge_idx[hei]%2 == 0){
        return m_rho_ve[hei->opposite()] + PI;
    }
    else{
        return m_rho_ve[hei->opposite()];
    }
}

double Connection::angleVertexToFacet(Halfedge_iterator hei){
    // head vertex
    return m_rho_vt[hei->next()];
}

Vector_3 Connection::rotate(Vector_3 v, Vector_3 axis, double rad){
    axis /= sqrt(axis.squared_length());
    Vector_3 vt = CGAL::cross_product(axis, v);
    
    return std::cos(rad)*v+std::sin(rad)*vt;
}

Vector_3 Connection::hevector(Halfedge_iterator hei){
    return hei->vertex()->point() - hei->opposite()->vertex()->point();
}

bool Connection::isEdge(Halfedge_iterator hei){
    if(m_halfedge_idx[hei]%2 == 0){
        return true;
    }
    else{
        return false;
    }
}

int Connection::edgeIdxFromHalfedge(Halfedge_iterator hei){
    if(m_halfedge_idx[hei]%2 == 0){
        return m_edge_idx[hei];
    }
    else{
        return m_edge_idx[hei->opposite()];
    }
}

Edge_iterator Connection::edgeFromHalfedge(Halfedge_iterator hei){
    if(m_halfedge_idx[hei]%2 == 0){
        return hei;
    }
    else{
        return hei->opposite();
    }
}

double Connection::cotan(Halfedge_iterator hei){
    Vector_3 v1 = hei->vertex()->point() - hei->next()->vertex()->point();
    Vector_3 v2 = hei->prev()->vertex()->point() - hei->next()->vertex()->point();
    
    double cos_theta = v1*v2/(norm(v1)*norm(v2));
    double sin_theta = norm(cross(v1,v2))/(norm(v1)*norm(v2));

    // hack!!!!!! tan is actually computed!!!!! change later ...
    return sin_theta / cos_theta;
}

double Connection::tan(Halfedge_iterator hei){
    Vector_3 v1 = hei->vertex()->point() - hei->next()->vertex()->point();
    Vector_3 v2 = hei->prev()->vertex()->point() - hei->next()->vertex()->point();
    
    double cos_theta = v1*v2/(norm(v1)*norm(v2));
    double sin_theta = norm(cross(v1,v2))/(norm(v1)*norm(v2));

    return sin_theta / cos_theta;
}

void Connection::frame(Vector_3& f1, Vector_3& f2, Vertex_iterator vi){
    f1 = hevector(vi->halfedge()->opposite()) - (hevector(vi->halfedge()->opposite())*normal(vi))*normal(vi);
    f1 /= sqrt(f1.squared_length());

    //f2 = CGAL::cross_product(normal(vi), f1);
    f2 = rotate(f1, normal(vi), PI/2);
}

Vector_3 Connection::vectorOnHalfedge(Halfedge_iterator hei, double a){
    Vector_3 v = hevector(hei);
    v /= norm(v);
    Vector_3 n = normal(hei);
    Vector_3 rv = rotate(v, n, a);
    return rv;
}

Vector_3 Connection::vectorOnFacet(Facet_iterator fi, double a){
    Vector_3 v = hevector(fi->halfedge());
    v /= norm(v);
    Vector_3 n = normal(fi);
    Vector_3 rv = rotate(v, n, a);
    return rv;
}

/******************************************************\
 * 
 * Debug Fuctions
 * 
\******************************************************/

void Connection::computeOptimalHalfedgeBasedConnection(){
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        double lcc = m_levi_civita[hei];
        double c1 = m_rho_vt[hei] - m_rho_ve[hei];
        double c2 = m_rho_vt[hei->next()] - m_rho_ve[hei->opposite()] - PI - 2*PI*m_n[hei->opposite()];
        //std::cout<<lcc<<" "<<c1<<" "<<c2<<std::endl;
        m_oc[hei] = (c1+c2)/2;
    }
}

void Connection::testBorderEdge(){
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        std::cout<<"Is edge representative border halfedge? "<<ei->is_border()<<std::endl;
        if(ei->is_border()){
            std::cout<<"WTF! edge representative maybe border halfedge ..."<<std::endl;
        }
    }
}

void Connection::testLeviCivitaConnection(){
    // intialize spannign tree
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        m_dinst[ei] = false;
    }

    // breath first search
    std::map<Facet_iterator, bool> visited;
    std::map<Facet_iterator, bool> placed;
    for(Facet_iterator fi = m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        visited[fi] = false;
        placed[fi] = false;
    }

    std::list<Facet_iterator> toBeVisited;
    toBeVisited.push_back(m_mesh.facets_begin());
    placed[m_mesh.facets_begin()] = true;
    m_lc_field[m_mesh.facets_begin()] = 0;

    while(!toBeVisited.empty()){
        //std::cout<<toBeVisited.size()<<std::endl;
        Facet_iterator fi = toBeVisited.front();
        toBeVisited.pop_front();
        m_lc_vf[fi] = vectorOnFacet(fi, m_lc_field[fi]);
        visited[fi] = true;

        Halfedge_iterator hei = fi->halfedge();
        int k = 0;
        do{
            if(!visited[hei->opposite()->facet()] && !placed[hei->opposite()->facet()]){
                ++k;
                toBeVisited.push_back(hei->opposite()->facet());
                placed[hei->opposite()->facet()] = true;

                Edge_iterator ei = hei;
                if(!isEdge(ei)){
                    ei = ei->opposite();
                }
                m_dinst[ei] = true;

                //m_lc_field[hei->opposite()->facet()] = m_lc_field[fi] + m_levi_civita[hei] - m_levi_civita[hei->opposite()] + PI;
                m_lc_field[hei->opposite()->facet()] = m_lc_field[fi] + m_oc[hei] - m_oc[hei->opposite()] + PI;
            }
            hei = hei->next();
        }while(hei != fi->halfedge());
        //std::cout<<k<<std::endl;
    }

    // write
    std::ofstream out("f_lc_frames.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Frames per simplex"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    int numc = m_mesh.size_of_facets();

    out<<"POINTS "<<2*numc<<" double"<<std::endl;
    for(Facet_iterator fi = m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        Point_3 p1 = baryCenter(fi);
        Vector_3 hev = m_lc_vf[fi];
        Point_3 p2 = p1 + 0.2*m_ael*hev/norm(hev);
        out<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<std::endl;
        out<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<std::endl;
    }

    out<<"CELLS "<<2*numc<<" "<<5*numc<<std::endl;
    for(int i=0; i<numc; ++i){
        out<<"2 "<<2*i<<" "<<2*i+1<<std::endl;
    }
    for(int i=0; i<numc; ++i){
        out<<"1 "<<2*i<<std::endl;
    }

    out<<"CELL_TYPES "<<2*numc<<std::endl;
    for(int i=0; i<numc; ++i){
        out<<"3"<<std::endl;
    }
    for(int i=0; i<numc; ++i){
        out<<"1"<<std::endl;
    }

    out.close();
}

void Connection::testInnerProduct(){
    std::cout<<"1 form same: "<<std::endl;
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        std::cout<<ipSame1F(hei)<<" "<<area(hei->facet())<<std::endl;
    }
    std::cout<<"1 form diff: "<<std::endl;
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        std::cout<<ipDiff1F(hei)<<std::endl;
    }

    std::cout<<"0 form same: "<<std::endl;
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        std::cout<<ipSame0F(hei)<<std::endl;
    }
    std::cout<<"0 form diff: "<<std::endl;
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        std::cout<<ipDiff0F(hei)<<std::endl;
    }
}

void Connection::testCotan(){
    std::cout<<"Testing cotan ..."<<std::endl;
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        std::cout<<cotan(hei)<<std::endl;
    }
}

void Connection::testAngleDefect(){
    std::cout<<"Testing angle defect"<<std::endl;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        std::cout<<defect(vi)<<std::endl;
    }
}

void Connection::testOrientedAngle(){
    std::cout<<"Testing oriented angle computation ..."<<std::endl;
    std::cout<<"Target: 1.047 (1/3*PI) # "<<"Result: "<<orientedAngle(Vector_3(1,0,0), Vector_3(1,1.732,0), Vector_3(0,0,1))<<std::endl;
    std::cout<<"Target: (2/3*PI) # "<<"Result: "<<orientedAngle(Vector_3(1,0,0), Vector_3(-1,1.732,0), Vector_3(0,0,1))<<std::endl;
    std::cout<<"Target: (-PI) # "<<"Result: "<<orientedAngle(Vector_3(1,0,0), Vector_3(-1,0,0), Vector_3(0,0,1))<<std::endl;
    std::cout<<"Target: (-2/3*PI) # "<<"Result: "<<orientedAngle(Vector_3(1,0,0), Vector_3(-1,-1.732,0), Vector_3(0,0,1))<<std::endl;
    std::cout<<"Target: (-1/3*PI) # "<<"Result: "<<orientedAngle(Vector_3(1,0,0), Vector_3(1,-1.732,0), Vector_3(0,0,1))<<std::endl;
}

void Connection::testHalfedgeRotation(){
    std::map<Halfedge_iterator, double> heiro;

    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        //if(valence(hei->opposite()->vertex()) == 6){
        //    continue;
        //}

        
        double angle_vieij, angle_e, angle_vjeij;
        double epsilon;
        if(isEdge(hei)){
            angle_vieij = m_rho_ve[hei];
            angle_e = m_rho_e[hei];
            angle_vjeij = (m_rho_ve[hei->opposite()] + PI + 2*PI*m_n[hei->opposite()]);
            epsilon = - angle_vieij + angle_e + angle_vjeij;
        }
        else{
            angle_vieij = m_rho_ve[hei];
            angle_e = - m_rho_e[hei->opposite()] + 2*PI*m_q[hei];
            angle_vjeij = (m_rho_ve[hei->opposite()] + PI + 2*PI*m_n[hei->opposite()]);
            epsilon = - angle_vieij + angle_e + angle_vjeij;
        }

        heiro[hei] = epsilon;
    }

    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        std::cout<<heiro[ei]<<" "<<heiro[ei->opposite()]<<std::endl;
    }
}

void Connection::testFacetRotation(){
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        double angle_e;
        if(isEdge(hei)){
            angle_e = m_rho_e[hei];
        }
        else{
            angle_e = - m_rho_e[hei->opposite()] + 2*PI*m_q[hei];
        }

        double tau = -m_rho_vt[hei] + angle_e + m_rho_vt[hei->next()];
        std::cout<<tau<<std::endl;
    }
}

void Connection::checkDT(){
    std::cout<<"Checking DT ..."<<std::endl;
    double tg = 0;
    for(Facet_iterator fi = m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        Halfedge_iterator hei = fi->halfedge();
        double rho_1 = m_rho_e[edgeFromHalfedge(hei)];
        if(!isEdge(hei)){
            rho_1 = -rho_1 + m_q[hei]*2*PI;
        }
        double rho_2 = m_rho_e[edgeFromHalfedge(hei->next())];
        if(!isEdge(hei->next())){
            rho_2 = -rho_2 + m_q[hei->next()]*2*PI;
        }
        double rho_3 = m_rho_e[edgeFromHalfedge(hei->prev())];
        if(!isEdge(hei->prev())){
            rho_3 = -rho_3 + m_q[hei->prev()]*2*PI;
        }

        double t1 = -m_rho_vt[hei] + rho_1 + m_rho_vt[hei->next()];
        double t2 = -m_rho_vt[hei->next()] + rho_2 + m_rho_vt[hei->prev()];
        double t3 = -m_rho_vt[hei->prev()] + rho_3 + m_rho_vt[hei];

        double t = t1+t2+t3;
        //std::cout<<t1+t2+t3<<std::endl;
        if(t>0){
            std::cout<<"WTF !!!!"<<std::endl;
        }
        tg += t;
    }
    std::cout<<"Total Gauss curvature: "<<tg<<std::endl;
}

void Connection::checkDE(){
    std::cout<<"Checking DE ..."<<std::endl;
    for(Halfedge_iterator hei = m_mesh.halfedges_begin(); hei != m_mesh.halfedges_end(); ++hei){
        if(valence(hei->opposite()->vertex()) == 6){
            continue;
        }

        double lcc = m_levi_civita[hei];
        double c1 = m_rho_vt[hei] - m_rho_ve[hei];
        double c2 = m_rho_vt[hei->next()] - m_rho_ve[hei->opposite()] - PI - 2*PI*m_n[hei->opposite()];
        std::cout<<lcc-c1<<std::endl;
        //std::cout<<lcc<<" "<<c1<<" "<<c2<<std::endl;
    }
}

void Connection::checkLCCInducedDefect(){
    std::cout<<"Checking Levi-Civita connection induced angle defect ..."<<std::endl;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        Halfedge_around_vertex_circulator hcir = vi->vertex_begin();
        double ad = 0;
        double jump = 0;
        do{
            ad += m_levi_civita[hcir];
            ad -= m_levi_civita[hcir->opposite()];
            
            jump += 0.5+m_n[hcir->opposite()];

            ++hcir;
        }while(hcir != vi->vertex_begin());
        std::cout<<ad+jump*2*PI<<" "<<jump<<std::endl;
    }
}

void Connection::checkDefect(){
    std::cout<<"Checking angle defect induced by globally optimal connection ..."<<std::endl;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        Halfedge_around_vertex_circulator hcir = vi->vertex_begin();
        double ad = 0;
        double jump = 0;
        do{
            ad += m_rho_vt[hcir->next()] - m_rho_ve[hcir->opposite()] - PI - 2*PI*m_n[hcir->opposite()];
            ad -= m_rho_vt[hcir->opposite()] - m_rho_ve[hcir->opposite()];
            jump += 0.5+m_n[hcir->opposite()];

            ++hcir;
        }while(hcir != vi->vertex_begin());
        std::cout<<ad+jump*2*PI<<std::endl;
    }
}

void Connection::isSymmetric(SparseMatrix sm){
    int num = 0;
    int maxn = 10;
    std::cout<<std::setprecision(15);
    for (int k = 0; k < sm.outerSize(); ++k){
        for (SparseMatrix::InnerIterator iter(sm, k); iter; ++iter){
            if(std::fabs(iter.value() - sm.coeff(iter.col(), iter.row()))>EPS){
                std::cout<<"Error: Matrix is not symmetric ... "<<std::endl;
                std::cout<<"["<<iter.row()<<" "<<iter.col()<<"]: "<<" "<<iter.value()<<std::endl;
                std::cout<<"["<<iter.col()<<" "<<iter.row()<<"]: "<<" "<<sm.coeff(iter.col(), iter.row())<<std::endl;
                ++num;
            }

            if(num > maxn){
                std::cout<<"Non-symmetric entries exceeds number "<<maxn<<" ..."<<std::endl;
                std::abort();
            }
        }
    }
}

void Connection::countNonzeroEntries(SparseMatrix sm){
    std::cout<<"# nonzero entries: "<<sm.nonZeros()<<std::endl;
}

void Connection::assignLocalVE(){
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        Halfedge_iterator hei = vi->halfedge();

        double da = 2*PI - defect(vi);
        double ca = 0;
        do{
            Vector_3 v1 = hei->opposite()->vertex()->point() - vi->point();
            Vector_3 v2 = hei->next()->vertex()->point() - vi->point();
            ca += angle(v1, v2);
            
            d_local_ve_angle[hei->next()] = -2*PI*ca/da;

            hei = hei->next()->opposite();
        }while(hei != vi->halfedge());
    }
}

void Connection::assignLocalVT(){
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        //std::cout<<"Vertex ..."<<std::endl;
        Halfedge_iterator hei = vi->halfedge();

        double da = 2*PI - defect(vi);
        double ca = 0;
        do{
            d_local_vt_angle[hei->opposite()] = 2*PI*ca/da + m_levi_civita[hei->opposite()];

            Vector_3 v1 = hei->opposite()->vertex()->point() - vi->point();
            Vector_3 v2 = hei->opposite()->next()->vertex()->point() - vi->point();
            ca += angle(v1, v2);
            //std::cout<<da<<" "<<ca<<std::endl;
            hei = hei->opposite()->next()->next();
        }while(hei != vi->halfedge());
    }
}

void Connection::totalAngleDefect(){
    double td = 0;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        td += defect(vi);
    }
    std::cout<<"Total angle defect: "<<td<<std::endl;
}

void Connection::saveEnergyMatrix(){
    std::ofstream out("EnergyMatrix.txt");
    for(int i=0; i<m_A.outerSize(); ++i){
        for(SparseMatrix::InnerIterator iter(m_A,i); iter; ++iter){
            out<<iter.row()+1<<" "<<iter.col()+1<<" "<<iter.value()<<std::endl;
        }
    }
    out.close();
}

void Connection::saveVector(VectorXd v){
    std::ofstream out("vector.txt");
    for(int i=0; i<v.size(); ++i){
        out<<v(i)<<std::endl;
    }
    out.close();
} 

void Connection::saveFrames(){
    std::ofstream out("frames.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Frames per simplex"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    double numc = m_mesh.size_of_vertices() + m_mesh.size_of_halfedges()/2 + m_mesh.size_of_facets();

    out<<"POINTS "<<2*numc<<" double"<<std::endl;
    // vertex frame
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        Point_3 p1 = vi->point();
        Vector_3 hev = hevector(vi->halfedge()->opposite());
        Point_3 p2 = p1 + 0.2*m_ael*hev/norm(hev);
        out<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<std::endl;
        out<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<std::endl;
    }
    // edge frame
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        Point_3 p1 = ei->opposite()->vertex()->point() + (ei->vertex()->point() - ei->opposite()->vertex()->point())/2;
        Vector_3 hev = hevector(ei);
        Point_3 p2 = p1 + 0.2*m_ael*hev/norm(hev);
        out<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<std::endl;
        out<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<std::endl;
    }
    // facet frame
    for(Facet_iterator fi = m_mesh.facets_begin(); fi != m_mesh.facets_end(); ++fi){
        Point_3 p1 = baryCenter(fi);
        Vector_3 hev = hevector(fi->halfedge());
        Point_3 p2 = p1 + 0.2*m_ael*hev/norm(hev);
        out<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<std::endl;
        out<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<std::endl;
    }

    out<<"CELLS "<<numc<<" "<<3*numc<<std::endl;
    for(int i=0; i<numc; ++i){
        out<<"2 "<<2*i<<" "<<2*i+1<<std::endl;
    }

    out<<"CELL_TYPES "<<numc<<std::endl;
    for(int i=0; i<numc; ++i){
        out<<"3"<<std::endl;
    }

    out<<"CELL_DATA "<<numc<<std::endl;
    out<<"LOOKUP_TABLE C "<<numc<<std::endl;
    for(int i=0; i<m_mesh.size_of_vertices(); ++i){
        out<<"0.7 0.3 0.3 1"<<std::endl;
    }
    for(int i=0; i<m_mesh.size_of_halfedges()/2; ++i){
        out<<"0.3 0.7 0.3 1"<<std::endl;
    }
    for(int i=0; i<m_mesh.size_of_facets(); ++i){
        out<<"0.3 0.3 0.7 1"<<std::endl;
    }

    out.close();
}

void Connection::saveVertexToEdgeFrames(){
    assignLocalVE();

    std::ofstream out("ve_frames.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Frames per simplex"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    int numc = m_mesh.size_of_vertices() + m_mesh.size_of_halfedges();

    out<<"POINTS "<<2*numc<<" double"<<std::endl;
    // vertex frame
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        Point_3 p1 = vi->point();
        Vector_3 hev = hevector(vi->halfedge()->opposite());
        Point_3 p2 = p1 + 0.2*m_ael*hev/norm(hev);
        out<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<std::endl;
        out<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<std::endl;

        Halfedge_around_vertex_circulator hcir = vi->vertex_begin();
        do{
            Point_3 pp1 = pointOnHalfedge(hcir, 0.75);
            double a = m_rho_ve[hcir->opposite()];
            //double a = d_local_ve_angle[hcir->opposite()];
            Halfedge_iterator hei = hcir->opposite();
            if(m_halfedge_idx[hei]%2 != 0){
                hei = hei->opposite();
                a += PI;
            }
            Vector_3 vv = vectorOnHalfedge(hei, -a);
            vv /= norm(vv);
            Point_3 pp2 = pp1 + 0.25*m_ael*vv;

            out<<pp1.x()<<" "<<pp1.y()<<" "<<pp1.z()<<std::endl;
            out<<pp2.x()<<" "<<pp2.y()<<" "<<pp2.z()<<std::endl;

            ++hcir;
        }while(hcir != vi->vertex_begin());
    }

    out<<"CELLS "<<2*numc<<" "<<5*numc<<std::endl;
    for(int i=0; i<numc; ++i){
        out<<"2 "<<2*i<<" "<<2*i+1<<std::endl;
    }
    for(int i=0; i<numc; ++i){
        out<<"1 "<<2*i<<std::endl;
    }

    out<<"CELL_TYPES "<<2*numc<<std::endl;
    for(int i=0; i<numc; ++i){
        out<<"3"<<std::endl;
    }
    for(int i=0; i<numc; ++i){
        out<<"1"<<std::endl;
    }

    out.close();
}

void Connection::saveVertexToTriangleFrames(){
    assignLocalVT();

    std::ofstream out("vt_frames.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Frames per simplex"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    int numc = m_mesh.size_of_vertices() + 3*m_mesh.size_of_facets();

    out<<"POINTS "<<2*numc<<" double"<<std::endl;
    // vertex frame
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        Point_3 p1 = vi->point();
        Vector_3 hev = hevector(vi->halfedge()->opposite());
        Point_3 p2 = p1 + 0.2*m_ael*hev/norm(hev);
        out<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<std::endl;
        out<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<std::endl;

        Halfedge_around_vertex_circulator hcir = vi->vertex_begin();
        do{
            Point_3 pp1 = pointOnFacet(hcir->opposite(), 0.8);
            double a = m_rho_vt[hcir->opposite()];
            //double a = d_local_vt_angle[hcir->opposite()];
            Vector_3 vv = vectorOnFacet(hcir->opposite()->facet(), -a);
            vv /= norm(vv);
            Point_3 pp2 = pp1 + 0.25*m_ael*vv;

            out<<pp1.x()<<" "<<pp1.y()<<" "<<pp1.z()<<std::endl;
            out<<pp2.x()<<" "<<pp2.y()<<" "<<pp2.z()<<std::endl;

            ++hcir;
        }while(hcir != vi->vertex_begin());
    }

    out<<"CELLS "<<2*numc<<" "<<5*numc<<std::endl;
    for(int i=0; i<numc; ++i){
        out<<"2 "<<2*i<<" "<<2*i+1<<std::endl;
    }
    for(int i=0; i<numc; ++i){
        out<<"1 "<<2*i<<std::endl;
    }

    out<<"CELL_TYPES "<<2*numc<<std::endl;
    for(int i=0; i<numc; ++i){
        out<<"3"<<std::endl;
    }
    for(int i=0; i<numc; ++i){
        out<<"1"<<std::endl;
    }

    out.close();
}

void Connection::saveVertexSpanningTree(){
    std::ofstream out("vertex_st.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Frames per simplex"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    int nume = 0;
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        if(m_inst[ei]){
            ++nume;
        }
    }

    out<<"POINTS "<<2*nume<<" double"<<std::endl;
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        if(!m_inst[ei]){
            continue;
        }

        Point_3 p1 = ei->vertex()->point();
        Point_3 p2 = ei->opposite()->vertex()->point();
        out<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<std::endl;
        out<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<std::endl;
    }

    out<<"CELLS "<<nume<<" "<<3*nume<<std::endl;
    for(int i=0; i<nume; ++i){
        out<<"2 "<<2*i<<" "<<2*i+1<<std::endl;
    }

    out<<"CELL_TYPES "<<nume<<std::endl;
    for(int i=0; i<nume; ++i){
        out<<3<<std::endl;
    }
    out.close();
}

void Connection::saveFacetSpanningTree(){
    std::ofstream out("facet_st.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Facet spanning tree"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    int num = 0;
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        if(m_dinst[ei]){
            ++num;
        }
    }

    out<<"POINTS "<<3*num<<" double"<<std::endl;
    for(Edge_iterator ei = m_mesh.edges_begin(); ei != m_mesh.edges_end(); ++ei){
        if(!m_dinst[ei]){
            continue;
        }

        Point_3 p1 = baryCenter(ei->facet());
        Point_3 p2 = baryCenter(ei);
        Point_3 p3 = baryCenter(ei->opposite()->facet());
        out<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<std::endl;
        out<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<std::endl;
        out<<p3.x()<<" "<<p3.y()<<" "<<p3.z()<<std::endl;
    }

    out<<"CELLS "<<2*num<<" "<<3*2*num<<std::endl;
    for(int i=0; i<num; ++i){
        out<<"2 "<<3*i<<" "<<3*i+1<<std::endl;
        out<<"2 "<<3*i+1<<" "<<3*i+2<<std::endl;
    }

    out<<"CELL_TYPES "<<2*num<<std::endl;
    for(int i=0; i<num; ++i){
        out<<3<<std::endl;
        out<<3<<std::endl;
    }

    out.close();
}

void Connection::saveVertexOneRingRotation(){
    std::ofstream out("vertex_rotation.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Vertex onering rotation"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    int num = 0;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        //if(valence(vi) == 6){
        //    continue;
        //}

        if(m_vertex_idx[vi]%20 != 0 || valence(vi) != 6){
            continue;
        }

        num += valence(vi)+1;
    }

    out<<"POINTS "<<2*num<<" double"<<std::endl;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        //if(valence(vi) == 6){
        //    continue;
        //}

        if(m_vertex_idx[vi]%20 != 0 || valence(vi) != 6){
            continue;
        }

        Point_3 p1 = vi->point();
        Vector_3 f1, f2;
        frame(f1, f2, vi);
        Point_3 p2 = p1+0.25*m_ael*f1;
        out<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<std::endl;
        out<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<std::endl;

        Halfedge_around_vertex_circulator hcir = vi->vertex_begin();
        do{
            Point_3 pp1 = hcir->opposite()->vertex()->point();
            
            double ra;
            if(isEdge(hcir)){
                ra = - (-m_rho_e[hcir] + m_q[hcir->opposite()]*2*PI);
            }
            else{
                ra = -m_rho_e[hcir->opposite()];
            }

            Vector_3 ff1, ff2;
            frame(ff1, ff2, hcir->opposite()->vertex());
            Vector_3 vv = rotate(ff1, normal(hcir->opposite()->vertex()), ra);
            Point_3 pp2 = pp1+0.25*m_ael*vv;
            out<<pp1.x()<<" "<<pp1.y()<<" "<<pp1.z()<<std::endl;
            out<<pp2.x()<<" "<<pp2.y()<<" "<<pp2.z()<<std::endl;

            ++hcir;
        }while(hcir != vi->vertex_begin());
    }

    out<<"CELLS "<<num<<" "<<3*num<<std::endl;
    for(int i=0; i<num; ++i){
        out<<"2 "<<2*i<<" "<<2*i+1<<std::endl;
    }

    out<<"CELL_TYPES "<<num<<std::endl;
    for(int i=0; i<num; ++i){
        out<<3<<std::endl;
    }

    out.close();

}

void Connection::saveVertexOneRingUniversalTest(){
    std::ofstream out("one_ring.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Frames per simplex"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    //int numc = m_mesh.size_of_vertices() + m_mesh.size_of_halfedges();
    int num = 0;
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        if(valence(vi) == 6){
            continue;
        }

        num += valence(vi)+1;
    }

    out<<"POINTS "<<2*num<<" double"<<std::endl;
    // vertex frame
    for(Vertex_iterator vi = m_mesh.vertices_begin(); vi != m_mesh.vertices_end(); ++vi){
        if(valence(vi) == 6){
            continue;
        }

        Point_3 p1 = vi->point();
        Vector_3 hev = hevector(vi->halfedge()->opposite());
        Point_3 p2 = p1 + 0.2*m_ael*hev/norm(hev);
        out<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<std::endl;
        out<<p2.x()<<" "<<p2.y()<<" "<<p2.z()<<std::endl;

        Halfedge_around_vertex_circulator hcir = vi->vertex_begin();
        do{
            double a = m_rho_ve[hcir->opposite()];
            Halfedge_iterator hei = hcir->opposite();
            //if(m_halfedge_idx[hei]%2 != 0){
            //    hei = hei->opposite();
            //    a += PI;
            //}

            // apply epsilon rotation
            double epsilon;
            if(isEdge(hei)){
                epsilon = - m_rho_ve[hei] + m_rho_e[hei] + m_rho_ve[hei->opposite()] + PI + 2*PI*m_n[hei->opposite()];
            }
            else{
                epsilon = - m_rho_ve[hei]
                    + (- m_rho_e[hei->opposite()] + 2*PI*m_q[hei])
                    + m_rho_ve[hei->opposite()] + PI + 2*PI*m_n[hei->opposite()];
            }
            //std::cout<<"Epsilon: "<<epsilon<<std::endl;

            double aa = m_rho_ve[hei->opposite()] + PI + 2*PI*m_n[hei->opposite()];

            //Point_3 pp1 = pointOnHalfedge(hcir, 0.5);
            //Vector_3 vv = vectorOnHalfedge(hei, - a - 0.5*epsilon);
            //vv /= norm(vv);
            //Point_3 pp2 = pp1 + 0.25*m_ael*vv;

            Point_3 pp1 = hei->vertex()->point();
            Vector_3 ff1, ff2;
            frame(ff1, ff2, hei->vertex());
            Vector_3 vv = rotate(ff1, normal(hei->vertex()), - a - epsilon + aa);
            vv /= norm(vv);
            std::cout<<vv<<std::endl;
            Point_3 pp2 = pp1 + 0.25*m_ael*vv;

            out<<pp1.x()<<" "<<pp1.y()<<" "<<pp1.z()<<std::endl;
            out<<pp2.x()<<" "<<pp2.y()<<" "<<pp2.z()<<std::endl;

            ++hcir;
            //std::cout<<m_halfedge_idx[hcir]<<std::endl;
        }while(hcir != vi->vertex_begin());
    }

    out<<"CELLS "<<2*num<<" "<<5*num<<std::endl;
    for(int i=0; i<num; ++i){
        out<<"2 "<<2*i<<" "<<2*i+1<<std::endl;
    }
    for(int i=0; i<num; ++i){
        out<<"1 "<<2*i<<std::endl;
    }

    out<<"CELL_TYPES "<<2*num<<std::endl;
    for(int i=0; i<num; ++i){
        out<<"3"<<std::endl;
    }
    for(int i=0; i<num; ++i){
        out<<"1"<<std::endl;
    }

    out.close();
}