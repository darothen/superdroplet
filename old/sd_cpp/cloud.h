#ifdef CLOUD_H_
#define CLOUD_H_

class Cloud {
  public:
    
  private:
    std::vector<std::unique_ptr<Droplet>> droplets_;
}

#endif  // CLOUD_H_