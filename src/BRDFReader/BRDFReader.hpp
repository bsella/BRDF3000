#include "BRDFReader.h"

using namespace std;
using namespace Eigen;

namespace ChefDevr
{
	template <typename Scalar>
	BRDFReader<Scalar>::BRDFReader()
	{

	}

	template <typename Scalar>
	Matrix<Scalar> BRDFReader<Scalar>::create_brdfSet(const char *fileDirectory)
	{

	}

	template <typename Scalar>
	Vector<Scalar> BRDFReader<Scalar>::sample_brdf(const Vector<Scalar> &brdf, unsigned int num_sampling)
	{
		const unsigned int num_phi = 4 * num_sampling;
		const double stepTheta = 0.5 * M_PI / num_sampling;
		const double stepPhi = 2.0 * M_PI / num_phi;
		const unsigned int num_coefficients = num_sampling * num_phi * num_sampling * num_phi;
		Vector<Scalar> sampled_brdf{num_coefficients};
		unsigned int index = 0;

		for (unsigned int i = 0; i < num_sampling; i++){
			const double theta_in = i * stepTheta;

			for (unsigned int j = 0; j < num_phi; j++){
				const double phi_in = j * stepPhi;

				for (unsigned int k = 0; k < num_sampling; k++){
					const double theta_out = k * stepTheta;

					for (unsigned int l = 0; l < num_phi; l++){
						const double phi_out = l * stepPhi;
						lookup_brdf_val(brdf, theta_in, phi_in, theta_out, phi_out, sampled_brdf[index], sampled_brdf[index + 1], sampled_brdf[index + 2]);
						index += 3;
					}
				}
			}
		}

		return sampled_brdf;
	}

	template <typename Scalar>
	Vector<Scalar> BRDFReader<Scalar>::read_brdf(const char *filePath)
	{
		const FILE *file = fopen(filePath, "rb");
		if (!file) {
			throw BRDFReaderError{string{"The file "} + filePath + " could not have been opened"};
		}

		int dims[3];
		const unsigned int num_valuesNeeded = samplingResolution_thetaH *
								 samplingResolution_thetaD * samplingResolution_phiD / 2;

		fread(dims, sizeof(int), 3, file);
		const unsigned int num_values = dims[0] * dims[1] * dims[2];
		if (num_values != num_valuesNeeded){
				fclose(file);
				throw BRDFReaderError{string{"Dimensions don't match : "} + to_string(num_values) +
									  " is not equal to " + to_string(num_valuesNeeded)};
		}

		const unsigned int num_coefficents = 3 * num_values;
		Vector<Scalar> brdf{num_coefficents};
		fread(brdf.data(), sizeof(Scalar), num_coefficents, file);

		fclose(file);

		return brdf;
	}

	template <typename Scalar>
	void BRDFReader<Scalar>::lookup_brdf_val(const Vector<Scalar> &brdf, double theta_in, double phi_in,
						 					 double theta_out, double phi_out,
											 double& red_value,double& green_value,double& blue_value)
	{
		// Convert to halfangle / difference angle coordinates
		double theta_half, phi_half, theta_diff, phi_diff;

		std_coords_to_half_diff_coords(theta_in, phi_in, theta_out, phi_out,
									   theta_half, phi_half, theta_diff, phi_diff);


		// Find index.
		// Note that phi_half is ignored, since isotropic BRDFs are assumed
		const unsigned int index = phi_diff_index(phi_diff) +
				  theta_diff_index(theta_diff) * samplingResolution_phiD * 0.5 +
				  theta_half_index(theta_half) * samplingResolution_phiD * 0.5 * samplingResolution_thetaD;

		const double stepBlue = samplingResolution_thetaH * samplingResolution_thetaD * samplingResolution_phiD;

		red_value = brdf[index] * red_scale;
		green_value = brdf[index + stepBlue * 0.5] * green_scale;
		blue_value = brdf[index + stepBlue] * blue_scale;

		if (red_value < 0.0 || green_value < 0.0 || blue_value < 0.0) {
			throw BRDFReaderError{"Below horizon."};
		}
	}

    template <typename Scalar>
    void BRDFReader<Scalar>::std_coords_to_half_diff_coords(double theta_in, double phi_in, double theta_out, double phi_out,
                                        double& theta_half,double& phi_half,double& theta_diff,double& phi_diff)
    {
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


    template <typename Scalar>
    Vector3d BRDFReader<Scalar>::rotate_vector(const Vector3d &vector, const Vector3d &axis, double angle)
    {
        double cos_angle = cos(angle);
        Vector3d out = vector * cos_angle;

        double temp = axis.dot(vector);
        temp = temp * (1.0 - cos_angle);
        out += axis * temp;

        Vector3d cross = axis.cross(vector);
        out += cross * sin(angle);

        return out;
    }


    template <typename Scalar>
    Vector3d BRDFReader<Scalar>::compute_direction (double theta, double phi) {
        Vector3d direction;
        const double proj_vec = sin(theta);

        direction.x() = proj_vec * cos(phi);
        direction.y() = proj_vec * sin(phi);
        direction.z() = cos(theta);

        return direction.normalized();
    }


	template <typename Scalar>
	unsigned int BRDFReader<Scalar>::theta_half_index(double theta_half)
	{
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


	template <typename Scalar>
	unsigned int BRDFReader<Scalar>::theta_diff_index(double theta_diff)
	{
		int result = int(theta_diff / (M_PI * 0.5) * samplingResolution_thetaD);

		if (result < 0) {
			result = 0;
		} else if (result > samplingResolution_thetaD - 1) {
			result = samplingResolution_thetaD - 1;
		}

		return  (unsigned int)result;
	}


	template <typename Scalar>
	unsigned int BRDFReader<Scalar>::phi_diff_index(double phi_diff) {
		// Because of reciprocity, the BRDF is unchanged under
		// phi_diff -> phi_diff + M_PI
		if (phi_diff < 0.0) {
			phi_diff += M_PI;
		}

		int result = int(phi_diff / M_PI * samplingResolution_phiD / 2);

		if (result < 0) {
			result = 0;
		} else if (result > samplingResolution_phiD / 2 - 1) {
			result = samplingResolution_phiD / 2 - 1;
		}

		return (unsigned int)result;
	}


} // namespace ChefDevr