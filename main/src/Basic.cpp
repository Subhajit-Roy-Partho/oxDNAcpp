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

