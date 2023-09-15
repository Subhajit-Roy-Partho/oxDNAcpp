// Print out the full vector definition
template <typename S>
std::ostream &std::operator<<(std::ostream &os, const std::vector<S> &vector)
{
  for (auto element : vector)
  {
    os << element << " ";
  }
  return os;
}

class Particle{
public:
  int id, color, strand;
  LR_vector r = {0, 0, 0}, a1, a3;
  inline void operator&=(Particle S){
    r=S.r;
    a1=S.a1;
    a3=S.a3;
    // color=S.color;
    // strand=S.strand;
  }
  inline void operator+=(Particle S){
        r=S.r;
    a1=S.a1;
    a3=S.a3;
    color=S.color;
    strand=S.strand;
  }
};

struct Traj{
  LR_vector r;
  LR_vector a1;
  LR_vector a3;
};

template <typename A, std::size_t N>
A npMean(A (&vector)[N])
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