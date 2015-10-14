/**
 * @brief TSDF::rayCasting returns all the voxels that are traversed by a ray going from start to end
 * @param start : continous world position wheray_startre the ray starts
 * @param end   : continous world position where the ray end
 * @return vector of voxel ids hit by the ray in temporal order
 *
 * J. Amanatides, A. Woo. A Fast Voxel Traversal Algorithm for Ray Tracing. Eurographics '87
 */
 
std::vector<Eigen::Vector3i> TSDF::rayCasting(Eigen::Vector3d ray_start, Eigen::Vector3d ray_end) {
  std::vector<Eigen::Vector3i> visited_voxels;

  // This id of the first/current voxel hit by the ray.
  // Using floor (round down) is actually very important,
  // the implicit int-casting will round up for negative numbers.
  Eigen::Vector3i current_voxel(std::floor(ray_start[0]/_bin_size),
                                std::floor(ray_start[1]/_bin_size),
                                std::floor(ray_start[2]/_bin_size));

  // The id of the last voxel hit by the ray.
  // TODO: what happens if the end point is on a border?
  Eigen::Vector3i last_voxel(std::floor(ray_end[0]/_bin_size),
                             std::floor(ray_end[1]/_bin_size),
                             std::floor(ray_end[2]/_bin_size));

  // Compute normalized ray direction.
  Eigen::Vector3d ray = ray_end-ray_start;
  //ray.normalize();

  // In which direction the voxel ids are incremented.
  double stepX = (ray[0] >= 0) ? 1:-1; // correct
  double stepY = (ray[1] >= 0) ? 1:-1; // correct
  double stepZ = (ray[2] >= 0) ? 1:-1; // correct

  // Distance along the ray to the next voxel border from the current position (tMaxX, tMaxY, tMaxZ).
  double next_voxel_boundary_x = (current_voxel[0]+stepX)*_bin_size; // correct
  double next_voxel_boundary_y = (current_voxel[1]+stepY)*_bin_size; // correct
  double next_voxel_boundary_z = (current_voxel[2]+stepZ)*_bin_size; // correct

  // tMaxX, tMaxY, tMaxZ -- distance until next intersection with voxel-border
  // the value of t at which the ray crosses the first vertical voxel boundary
  double tMaxX = (ray[0]!=0) ? (next_voxel_boundary_x - ray_start[0])/ray[0] : DBL_MAX; //
  double tMaxY = (ray[1]!=0) ? (next_voxel_boundary_y - ray_start[1])/ray[1] : DBL_MAX; //
  double tMaxZ = (ray[2]!=0) ? (next_voxel_boundary_z - ray_start[2])/ray[2] : DBL_MAX; //

  // tDeltaX, tDeltaY, tDeltaZ --
  // how far along the ray we must move for the horizontal component to equal the width of a voxel
  // the direction in which we traverse the grid
  // can only be FLT_MAX if we never go in that direction
  double tDeltaX = (ray[0]!=0) ? _bin_size/ray[0]*stepX : DBL_MAX;
  double tDeltaY = (ray[1]!=0) ? _bin_size/ray[1]*stepY : DBL_MAX;
  double tDeltaZ = (ray[2]!=0) ? _bin_size/ray[2]*stepZ : DBL_MAX;

  if (current_voxel[0]!=last_voxel[0] && ray[0]<0) {current_voxel[0]--;}
  if (current_voxel[1]!=last_voxel[1] && ray[1]<0) {current_voxel[1]--;}
  if (current_voxel[2]!=last_voxel[2] && ray[2]<0) {current_voxel[2]--;}

  //Eigen::Vector3i first_voxel = current_voxel;
  visited_voxels.push_back(current_voxel);

  while(last_voxel != current_voxel) {
    if (tMaxX < tMaxY) {
      if (tMaxX < tMaxZ) {
        current_voxel[0] += stepX;
        tMaxX += tDeltaX;
      } else {
        current_voxel[2] += stepZ;
        tMaxZ += tDeltaZ;
      }
    } else {
      if (tMaxY < tMaxZ) {
        current_voxel[1] += stepY;
        tMaxY += tDeltaY;
      } else {
        current_voxel[2] += stepZ;
        tMaxZ += tDeltaZ;
      }
    }

    visited_voxels.push_back(current_voxel);
    /***
     * For debugging: if the voxels along the ray become too many,
     * we probably passed the targeted measurement.
     *
    if (visited_voxels.size() == 5000) {
      std::cout << std::endl;
      std::cout << "TOO MANY ELEMENTS! STOPING" << std::endl;
      std::cout << "First voxel: " << first_voxel[0] << "," << first_voxel[1] << "," << first_voxel[2] << std::endl;
      std::cout << "Last voxel: " << last_voxel[0] << "," << last_voxel[1] << "," << last_voxel[2] << std::endl;
      std::cout << "Ray: " << ray[0] << ", " << ray[1] << ", " << ray[2] << std::endl;
      std::cout << "Ray Start: " << ray_start[0] << ", " << ray_start[1] << ", " << ray_start[2] << std::endl;
      std::cout << "Ray End: " << ray_end[0] << ", " << ray_end[1] << ", " << ray_end[2] << std::endl;
      break;
    }*/
  }
  return visited_voxels;
}
