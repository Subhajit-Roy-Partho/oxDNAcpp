class Forces
{
public:
  std::vector<int> particles;
  std::vector<float> stiff, rate, position;
  std::vector<std::string> type;
  std::vector<LR_vector> dir, pos0;

  bool add(string name, int particle = -1, float stiff = 0, LR_vector normal = {0, 0, 1}, LR_vector pos = {0, 0, 0})
  {
    type.push_back(name);
    particles.push_back(particle);
    dir.push_back(normal);
    pos0.push_back(pos);
    this->stiff.push_back(stiff);
    return true;
  }
  bool addRepulsion(int particles = -1, float stiff = 1, LR_vector dir = {0, 0, 1}, float position = 0)
  {
    type.push_back("repulsion_plane");
    this->particles.push_back(particles);
    this->stiff.push_back(stiff);
    this->dir.push_back(dir);
    this->position.push_back(position);
    return true;
  }
  bool addHarmonic(int particles = 0, float stiff = 1, float rate = 0, LR_vector pos0 = {0, 0, 0}, LR_vector dir = {1, 0, 0})
  {
    type.push_back("trap");
    this->particles.push_back(particles);
    this->stiff.push_back(stiff);
    this->rate.push_back(rate);
    this->pos0.push_back(pos0);
    this->dir.push_back(dir);
    return true;
  }

  bool addHarmonicFromParticle(Particle *particle, float stiff = 1, float rate = 0)
  {
    return addHarmonic(particle->id, stiff, rate, particle->r);
  }

  bool save(string filename = "external_forces.txt")
  {
    ofstream file(filename);
    if (!file.is_open())
      return false;
    for (int i = 0; i < type.size(); i++)
    {
      file << "{" << std::endl;
      if (type[i] == "trap")
      {
        file << "type = trap\n";
        file << "particle = " << particles[0] << endl;
        particles.erase(particles.begin());
        file << "pos0 = " << pos0[0].x<<","<<pos0[0].y<<","<<pos0[0].z << endl;
        pos0.erase(pos0.begin());
        file << "stiff = " << stiff[0] << endl;
        stiff.erase(stiff.begin());
        file << "rate = " << rate[0] << endl;
        rate.erase(rate.begin());
        file << "dir = " << dir[0].x<<","<<dir[0].y<<","<<dir[0].z << endl;
        dir.erase(dir.begin());
      }
      else if (type[i] == "repulsion_plane")
      {
        file << "type = repulsion_plane" << endl;
        file << "particle = " << particles[0] << endl;
        particles.erase(particles.begin());
        file << "stiff = " << stiff[0] << endl;
        stiff.erase(stiff.begin());
        file << "dir = " << dir[0].x<<","<<dir[0].y<<","<<dir[0].z << endl;
        dir.erase(dir.begin());
        file << "position = " << position[0] << endl;
        position.erase(position.begin());
      }
      file << "}\n";
    }
    file.close();
    return true;
  }
};