#include "BRDFReader.h"

namespace ChefDevr {

    BRDFReader::BRDFReader() {

    }

    void BRDFReader::std_coords_to_half_diff_coords(double theta_in, double phi_in, double theta_out, double phi_out, double& theta_half, double& phi_half, double& theta_diff, double& phi_diff) {
        const Vector3d in = compute_direction(theta_in, phi_in);
        const Vector3d out = compute_direction(theta_out, phi_out);

        Vector3d halfway = (in + out) * 0.5f;
        halfway.normalize();

        theta_half = acos(halfway[2]);
        phi_half = atan2(halfway[1], halfway[0]);

        const Vector3d bi_normal = Vector3d::UnitY();
        const Vector3d normal = Vector3d::UnitZ();

        // compute diff vector
        const Vector3d temp = rotate_vector(in, normal , -phi_half);
        const Vector3d diff = rotate_vector(temp, bi_normal, -theta_half);

        theta_diff = acos(diff[2]);
        phi_diff = atan2(diff[1], diff[0]);
    }

    Vector3d BRDFReader::rotate_vector(const Vector3d &vector, const Vector3d &axis, double angle) {
        double cos_angle = cos(angle);
        Vector3d out = vector * cos_angle;

        double temp = axis.dot(vector);
        temp = temp * (1.0 - cos_angle);
        out += axis * temp;

        Vector3d cross = axis.cross(vector);
        out += cross * sin(angle);

        return out;
    }

    Vector3d BRDFReader::compute_direction (double theta, double phi) {
        Vector3d direction;
        const double proj_vec = sin(theta);

        direction.x() = proj_vec * cos(phi);
        direction.y() = proj_vec * sin(phi);
        direction.z() = cos(theta);

        return direction.normalized();
    }

    unsigned int BRDFReader::theta_half_index(double theta_half) {
        int result;

        if (theta_half <= 0.0) {
            result = 0;
        } else {
            const double theta_half_deg = ((theta_half / (M_PI / 2.0)) * samplingResolution_thetaH);
            const double temp = sqrt(theta_half_deg * samplingResolution_thetaH);

            result = (int)temp;
            if (result < 0) {
                result = 0;
            }

            if (result >= samplingResolution_thetaH) {
                result = samplingResolution_thetaH - 1;
            }
        }
        return (unsigned int)result;
    }

} // namespace ChefDevr
