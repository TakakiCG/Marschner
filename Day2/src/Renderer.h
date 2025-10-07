//
// Created by kango on 2023/04/03.
//

#ifndef DAY_3_RENDERER_H
#define DAY_3_RENDERER_H


#include <vector>
#include "Body.h"
#include "Camera.h"
#include <random>
#include "Eigen/Dense"

class Renderer {
public:
    std::vector<Body> bodies;

    Camera camera;
    Color bgColor;

    /// 乱数生成器
    mutable std::mt19937_64 engine;
    mutable std::uniform_real_distribution<> dist;

    Renderer(const std::vector<Body> &bodies, Camera camera, Color bgColor=Color::Zero());

    double rand() const;

    bool hitScene(const Ray &ray, RayHit &hit) const;

    Image passTracingRender(const unsigned int &samples) const;

    Image passTracingRender_antiAreas(const unsigned int &samples) const;

    Color trace(const Ray &ray, const RayHit &hit) const;

    Color traceMarschner(const Ray &ray, const RayHit &hit) const;

    void diffuseSample(const Eigen::Vector3d &incidentPoint, const Eigen::Vector3d &normal, Ray &out_Ray) const;

    static void marschnerSampleHair(const Ray &in_Ray, const Eigen::Vector3d &incidentPoint, const Eigen::Vector3d &axis,
                             Ray &out_ray) ;

    double FresnelFunction(double cos_gammaI, double eta) const;

    static void computeLocalFrame(const Eigen::Vector3d &w, Eigen::Vector3d &u, Eigen::Vector3d &v);

    void debugSampleHFromCamera(int n) const;


    /// Visualize theta and phi
    struct RaySampleData {
        int    sampleIdx     = 0;    // 何番目のサンプルか
        double incidentTheta = 0.0;  // θ_in  [rad]
        double incidentPhi   = 0.0;
        double outgoingTheta = 0.0;  // θ_out [rad]
        double outgoingPhi   = 0.0;  // φ_out [rad]

        std::string hitType      = "Unknown";   // "Sphere" or "Cylinder"
        std::string scatterType  = "None";      // "Diffuse" or "Specular" (or "None")
        std::string reason      = "N/A";       // why scatterType == None
    };

    mutable std::vector<RaySampleData> g_samples;
    mutable std::mutex                 g_mutex;

    Color trace(const Ray &ray, const RayHit &hit, bool recordAngles, int sampleIdx) const;

    void diffuseSample(const Eigen::Vector3d &incidentPoint, const Eigen::Vector3d &normal, Ray &out_Ray, bool recordAngles, RaySampleData* rec) const;
    static void marschnerSampleHair(const Ray& in_Ray,
                                    const Eigen::Vector3d& incidentPoint,
                                    const Eigen::Vector3d& axis,
                                    Ray& out_Ray,
                                    bool recordAngles = false,
                                    RaySampleData* rec = nullptr);


/// Marschner
//    struct BSDFParams{
//        const double eta = 1.55;
//        const double alpha_R = 0;
//        const double beta_R = 0;
//    };

//    Color fur_bsdf(const Eigen::Vector3d &wi, const Eigen::Vector3d &wo, const BSDFParams &params) const;

//    void marschnerSample(const Ray &in_ray, const Eigen::Vector3d &incidentPoint,
//                         const Eigen::Vector3d &hairDir, const Material &material,
//                         Ray &out_Ray) const;
//
//    Color marschnerShading(const Eigen::Vector3d &V, const Eigen::Vector3d &L,
//                           const Eigen::Vector3d &H, const Material &material) const;

};


#endif //DAY_3_RENDERER_H
