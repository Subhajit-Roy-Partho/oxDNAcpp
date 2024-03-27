// Print out the full vector definition
template <typename S>
std::ostream &operator<<(std::ostream &os, const std::vector<S> &vector){
  for (auto element : vector)
  {
    os << element << " ";
  }
  return os;
}

template <typename T>std::string to_string_with_precision(const T a_value, const int n = 6){
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}

struct Patch {
    LR_vector position; //the position of the patch with respect to the CM of the particle
    LR_vector a1;  //vector that is to be compared against the vector connecting the patches r_pp, needs to be parallel
    LR_vector a2; // vector which needs to be parallel with a2 on the complementary patch
    int id=0; //the id of the patch; it is used during initialization to assign patches to particles according to input file; sets the type of patch
    int index ; //this is the unique index of the patch in the simulation
    bool active = false; //is the patch on or not
    int locked_to_particle=-1; //id of the particle this patch is bound to; used to make sure 1 patch is only in 1 interaction
    int locked_to_patch=-1; //id of the particle this patch is bound to; used to make sure 1 patch is only in 1 interaction
    number locked_energy=0;
    bool spherical = false; //is the patch in spherical coordinates or not

    int color=-1; //this is the color of the patch
    number strength=1;  //sets the strength of the interaction/
    number radius=0.4;
    // number a1_x=1, a1_y=0, a1_z=0;
    // number a2_x=1, a2_y=0, a2_z=0;
};

class Particle{
public:
  int id, color, strand=0,type=0;
  double radius;
  std::string name;
  LR_vector r = {0, 0, 0}, a1={1,1,0}, a3={0,0,-1};
  std::vector<int> connector; // ids of particles connected
  std::vector<double> spring,eqRadius; // saves the spring constant of the ids
  std::vector<int> patches;
  // inline void operator&=(Particle S){
  //   r=S.r;
  //   a1=S.a1;
  //   a3=S.a3;
  //   // color=S.color;
  //   // strand=S.strand;
  // }
  inline void operator+=(Particle S){
    r=S.r;
    a1=S.a1;
    a3=S.a3;
    color=S.color;
    strand=S.strand;
    connector=S.connector;
    spring=S.spring;
    eqRadius=S.eqRadius;
  }
};

struct Traj{
  unsigned int time;
  std::vector<LR_vector> r;
  std::vector<LR_vector> a1;
  std::vector<LR_vector> a3;
  void updateParticleNumber(int particleNum){
    r.resize(particleNum);
    r.shrink_to_fit();
    a1.resize(particleNum);
    a1.shrink_to_fit();
    a3.resize(particleNum);
    a3.shrink_to_fit();
  }
};

template <typename A, std::size_t N>A npMean(A (&vector)[N])
{
  A result = vector[0];
  for (int i = 1; i < N; i++)
  {
    result += vector[i];
  }
  return result / N;
};

template <typename A> double npNorm(A vector, int norm = 2){
  if (std::is_same<A, LR_vector>::value)
    return pow(pow(vector.x, norm) + pow(vector.y, norm) + pow(vector.z, norm), 1.0 / (double)norm);

  return 0;
}

template <typename A> A npFloor(A vector){
  if (std::is_same<A, LR_vector>::value)
    return (LR_vector){std::floor(vector.x), std::floor(vector.y), std::floor(vector.z)};
}

template <typename A>A npRound(A vector){
  if (std::is_same<A, LR_vector>::value)
    return (LR_vector){std::round(vector.x), std::round(vector.y), std::round(vector.z)};
}

// This is basic sgn function whre 1 for >0 and -1 for <0 and 0 for 0.
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}