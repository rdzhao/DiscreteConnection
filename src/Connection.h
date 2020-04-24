#ifndef __CONNECTION__
#define __CONNECTION__

#include "Utility.h"

/******************************************************\
 * 
 * Class: Connection
 * Dealing with discrete connection
 * 
 * Frame conventions: 
 * Vertex: along halfedge()->opposite()
 *  (Warning: vertex->halfedge() is the halfedge that points to this vertex !!!)
 * Edge: along the cooresponding halfedge
 * Facet: along the halfedge() 
 * 
\******************************************************/

class Connection
{
public:
    Connection();
    ~Connection();

    // reference accessor
    Polyhedron& Mesh();
    int VertexIdx(Vertex_iterator vi);
    int EdgeIdx(Edge_iterator ei);
    int HalfedgeIdx(Halfedge_iterator hei);
    int FacetIdx(Facet_iterator fi);
    double RhoE(Halfedge_iterator hei);
    double RhoVE(Halfedge_iterator hei);
    double RhoVT(Halfedge_iterator hei);

    // major functions
    void read(std::string mesh_fn);
    void preprocessing();
    void compute();
    void write();


protected:
    // major internal functions
    void readMesh(std::string mesh_fn);

    void indexElements();
    void computeAverageEdgeLength();

    void assembleLinearSystem();
    void solveLinearSystem();
    void computeFieldByConnection();
    void computeVertexBasedVF();

    void writeMeshToOFF();
    void writeMeshToVTK();
    void writeMeshToVTKWithField();

    // minor helper functions
    void computeLeviCivitaConnection();
    void determineIntegerJump();
    void assembleEnergyMatrix();
    void assembleEnergyVector();

    void determine_n(); // (high confidence)
    void determine_q(); // (high confidence)

    // Diagonal block: no symmetric filling,
    // because we fill line by line for each target.
    // Off-diagonal block: symmetric filling, because we only start from one type
    // Indexing: [v->t, v->e, v->v]
    void assembleRhovtRhovt_Et(std::vector<Triplet>& triplets); // high confidence
    void assembleRhovRhov_Et(std::vector<Triplet>& triplets); // high confidence
    void assembleRhovtRhov_Et(std::vector<Triplet>& triplets); // high confidence
    void assembleRhovtRhovt_Ee(std::vector<Triplet>& triplets); // lower confidence
    void assembleRhoveRhove_Ee(std::vector<Triplet>& triplets); // lower confidence
    void assembleRhovtRhove_Ee(std::vector<Triplet>& triplets); // lower confidence

    void assembleVectorFromEt(VectorXd& v); // lower confidence
    void assembleVectorFromEe(VectorXd& v); // lower confidence

    void fixGauge();

    double ipSame0F(Halfedge_iterator hei); // inner product int phi_i*phi_i on edge (high confidence)
    double ipDiff0F(Halfedge_iterator hei); // inner product int phi_i*phi_j on edge (high confidence)
    double ipSame1F(Halfedge_iterator hei); // inner product int phi_ij dot phi_ij on triangle (high confidence)
    double ipDiff1F(Halfedge_iterator hei); // inner product int phi_ij dot phi_jk on traingle (high confidence)

    // basic geometric properties
    double norm(Vector_3 v);
    Vector_3 cross(Vector_3 v1, Vector_3 v2);

    Vector_3 normal(Vertex_iterator vi);
    double defect(Vertex_iterator vi);
    int valence(Vertex_iterator vi);
    Point_3 baryCenter(Edge_iterator ei);
    Vector_3 normal(Edge_iterator ei);
    Vector_3 height(Halfedge_iterator hei);
    double internalAngle(Halfedge_iterator hei);
    double length(Halfedge_iterator hei);
    Vector_3 hevector(Halfedge_iterator hei); // halfedge vector
    double tan(Halfedge_iterator hei);
    double cotan(Halfedge_iterator hei);
    double area(Facet_iterator fi);
    Point_3 baryCenter(Facet_iterator fi);
    Vector_3 normal(Facet_iterator fi);

    // basic operations
    double angle(Vector_3 v1, Vector_3 v2);
    double orientedAngle(Vector_3 v1, Vector_3 v2, Vector_3 n);

    Point_3 pointOnHalfedge(Halfedge_iterator hei, double t); 
    Point_3 pointOnFacet(Halfedge_iterator hei, double t);

    double angleVertexToEdge(Halfedge_iterator hei);
    double angleVertexToFacet(Halfedge_iterator hei);

    Vector_3 rotate(Vector_3 v, Vector_3 axis, double rad);

    void frame(Vector_3& f1, Vector_3& f2, Vertex_iterator vi);
    Vector_3 vectorOnHalfedge(Halfedge_iterator hei, double a);
    Vector_3 vectorOnFacet(Facet_iterator fi, double a);

    bool isEdge(Halfedge_iterator hei);
    int edgeIdxFromHalfedge(Halfedge_iterator hei);
    Edge_iterator edgeFromHalfedge(Halfedge_iterator hei);
 
    // debug functions
    void computeOptimalHalfedgeBasedConnection();

    void testBorderEdge();
    void testLeviCivitaConnection();
    void testInnerProduct();
    void testCotan();
    void testAngleDefect();
    void testOrientedAngle();
    void testHalfedgeRotation();
    void testFacetRotation();

    void checkDT();
    void checkDE();
    void checkLCCInducedDefect();
    void checkDefect();

    void isSymmetric(SparseMatrix sm);
    void countNonzeroEntries(SparseMatrix sm);

    void assignLocalVE();
    void assignLocalVT();
    void totalAngleDefect();

    void saveEnergyMatrix();
    void saveVector(VectorXd v);
    void saveFrames();
    void saveVertexToEdgeFrames();
    void saveVertexToTriangleFrames();
    void saveVertexSpanningTree();
    void saveFacetSpanningTree();
    void saveVertexOneRingRotation();
    void saveVertexOneRingUniversalTest();

protected:
    Polyhedron m_mesh;

    double m_ael; // average edge length

    std::map<Vertex_iterator, int> m_vertex_idx;
    std::map<Edge_iterator, int> m_edge_idx;
    std::map<Halfedge_iterator, int> m_halfedge_idx;
    std::map<Facet_iterator, int> m_facet_idx;

    // discrete connection
    // CGAL's halfedge->vertex() points to the head vertex
    // but we prefer using tail to encode information
    // for v->t, we use tail vertex to point to the associated facet
    // for v->e, we use tail vertex to point to the associated halfedge
    std::map<Halfedge_iterator, double> m_levi_civita; // Levi-Civita connection
    std::map<Halfedge_iterator, double> m_rho_vt; // tail vertex to face
    std::map<Halfedge_iterator, double> m_rho_ve; // tail vertex to halfedge v_i to e_ij with orientation
    std::map<Edge_iterator, double> m_rho_e; // with edge orientation 

    // interger jumps stored on halfedges
    
    // base direction is from v_tail to e_{tail,head}
    // this integer for v_tail to e_{head,tail} 
    // refer to equation before section 4.2
    std::map<Halfedge_iterator, int> m_n;

    // base direction is edge direction
    // halfedge direction is \rho_ij direction
    // if halfedge is with edge direction, m_q is 0 (representing rho_ij)
    // if halfedge is flipped with edge direction, m_q is as follows (representing rho_ji)
    // \rho_flipped = -\rho + m_q*2*PI
    std::map<Halfedge_iterator, int> m_q;

    SparseMatrix m_SMEt;
    SparseMatrix m_SMEe;
    SparseMatrix m_A;
    VectorXd m_b;

    // generated fields
    std::map<Vertex_iterator, double> m_field;
    std::map<Vertex_iterator, Vector_3> m_vf;
    std::map<Facet_iterator, double> m_lc_field;
    std::map<Facet_iterator, Vector_3> m_lc_vf; // levi-civita vector field

    // data for debug
    std::map<Halfedge_iterator, double> d_local_ve_angle;
    std::map<Halfedge_iterator, double> d_local_vt_angle;

    // optimal halfedge based connection
    std::map<Halfedge_iterator, double> m_oc;

    // spanning tree for parallel transport along edge
    std::map<Edge_iterator, bool> m_inst; // in spanning tree
    std::map<Edge_iterator, bool> m_dinst; // dual in spanning tree for edge rotation angle 
};

#endif


