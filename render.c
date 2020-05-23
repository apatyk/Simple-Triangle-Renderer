//
//  Adam Patyk
//  render.c
//  ECE 6680 Lab 5: Triangle Rendering
//
//  Copyright © 2020 Adam Patyk. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define ROWS 256
#define COLS 256

typedef struct bounding_box_coords_type {
    float *left, *right, *top, *bottom, *topleft;
} bbox_coords_t;

void read_ply_header(FILE *, int *, int *);
void read_ply_data(FILE *, int, int, float **, int **);
float calc_bbox(float **, int, float *, float *, float *);
void calc_cam_vec(float *, float *, int *, float *, float);
bbox_coords_t *calc_3D_coords(float *, float *, float *, float);
void draw_image(unsigned char *, float *, int, int, float **, int, int **,  int, bbox_coords_t *, float *);
void write_image(unsigned char *, int, int, char *);
void vec_mat_mult(float *, float[3][3]);
float dot_product(float *, float *, int);
float *cross_product(float *, float *);

int main(int argc, char *const *argv) {
    FILE    *fpt;
    int     i, num_vertices, num_faces, **faces;
    float   E, **vertices, *z_buf;
    float center[3] = {0};
    float camera[3] = {1, 0, 0};
    float up[3]     = {0, 0, 1};
    unsigned char *img;
    bbox_coords_t *coordinates;

    if (argc != 5) {
        fprintf(stderr, "Usage: ./render <file.ply> <camera_x> <camera_y> <camera_z>\n");
        exit(0);
    }

    // organize vector of camera rotations from command line
    int cam_rotate[3] = {atoi(argv[2]), atoi(argv[3]), atoi(argv[4])};
    // read ply header from file
    fpt = fopen(argv[1], "rt");
    read_ply_header(fpt, &num_vertices, &num_faces);
    vertices = (float **)calloc(num_vertices, sizeof(float *));

    for (i = 0; i < num_vertices; i++)
        vertices[i] = (float *)calloc(3, sizeof(float));

    faces = (int **)calloc(num_faces, sizeof(int *));

    for (i = 0; i < num_faces; i++)
        faces[i] = (int *)calloc(3, sizeof(int));

    // read vertices and faces from file
    read_ply_data(fpt, num_vertices, num_faces, vertices, faces);
    // calculate bounding box
    float min[3] = {vertices[0][0], vertices[0][1], vertices[0][2]};
    float max[3] = {vertices[0][0], vertices[0][1], vertices[0][2]};
    E = calc_bbox(vertices, num_vertices, min, max, center);
    // calculate camera position and orientation
    calc_cam_vec(camera, up, cam_rotate, center, E);
    // calculate 3D coordiantes bounding the image
    coordinates = calc_3D_coords(camera, up, center, E);
    // draw the image
    img = (unsigned char *)calloc(ROWS * COLS, sizeof(unsigned char));
    z_buf = (float *)calloc(ROWS * COLS, sizeof(float));

    // initialize z buffer to all high values
    for (i = 0; i < ROWS * COLS; i++)
        z_buf[i] = FLT_MAX;

    draw_image(img, z_buf, ROWS, COLS, vertices, num_vertices, faces, num_faces, coordinates, camera);
    write_image(img, ROWS, COLS, "render.ppm");
    printf("Rendering complete\n");
    // cleanup
    fclose(fpt);

    for (i = 0; i < num_faces; i++)
        free(faces[i]);

    for (i = 0; i < num_vertices; i++)
        free(vertices[i]);

    free(coordinates->left);
    free(coordinates->right);
    free(coordinates->top);
    free(coordinates->bottom);
    free(coordinates->topleft);
    free(faces);
    free(vertices);
    free(coordinates);
    free(img);
    free(z_buf);
    return 0;
}

void read_ply_header(FILE *fpt, int *num_vertices, int *num_faces) {
    int i, ret;
    char file_type[4];
    char tmp[256] = {0};
    // read in file type and excess tags
    ret = fscanf(fpt, "%s %*s %*s %*s %*f ", file_type);

    if (ret == 0 || strcmp(file_type, "ply")) {
        fprintf(stderr, "Must be .ply file!\n");
        exit(0);
    }

    // read in number of vertices
    fscanf(fpt, "%*s %*s %d", num_vertices);

    for (i = 0; i < 3; i++)
        fscanf(fpt, "%*s %*s %*s");

    // read in number of faces
    fscanf(fpt, "%*s %*s %d", num_faces);

    // read in rest of header
    while (strcmp(tmp, "end_header"))
        fscanf(fpt, "%s", tmp);

    printf("vertices: %d, faces: %d\n", *num_vertices, *num_faces);
}

void read_ply_data(FILE *fpt, int num_vertices, int num_faces, float **vertices, int **faces) {
    int i, j;

    for (i = 0; i < num_vertices; i++)
        for (j = 0; j < 3; j++)
            (void)fscanf(fpt, "%f", &vertices[i][j]);

    for (i = 0; i < num_faces; i++)
        (void)fscanf(fpt, "%*d %d %d %d", &faces[i][0], &faces[i][1], &faces[i][2]);
}

float calc_bbox(float **vertices, int num_vertices, float *min, float *max, float *center) {
    int i, j;
    float E;
    float total[3] = {0};

    // find min & max x, y, z
    for (i = 0; i < num_vertices; i++) {
        for (j = 0; j < 3; j++) {
            if (vertices[i][j] < min[j]) min[j] = vertices[i][j];

            if (vertices[i][j] > max[j]) max[j] = vertices[i][j];

            total[j] += vertices[i][j];
        }
    }

    E = max[0] - min[0];

    for (j = 0; j < 3; j++) {
        // find center x, y, z
        center[j] = total[j] / num_vertices;

        // find extent of bounding box
        if (max[j] - min[j] > E) E = max[j] - min[j];
    }

    printf("max: x=%f, y=%f, z=%f\n", max[0], max[1], max[2]);
    printf("min: x=%f, y=%f, z=%f\n", min[0], min[1], min[2]);
    printf("center: x=%f, y=%f, z=%f\n", center[0], center[1], center[2]);
    printf("extent: %f\n", E);
    return E;
}

void calc_cam_vec(float *camera, float *up, int *cam_rotate, float *center, float E) {
    float theta_x = cam_rotate[0] * M_PI / 180.0;
    float theta_y = cam_rotate[1] * M_PI / 180.0;
    float theta_z = cam_rotate[2] * M_PI / 180.0;
    // rotation matrices
    float Rx[3][3] = {
        {1, 0, 0},
        {0, cos(theta_x), -sin(theta_x)},
        {0, sin(theta_x), cos(theta_x)}
    };
    float Ry[3][3] = {
        {cos(theta_y), 0, sin(theta_y)},
        {0, 1, 0},
        { -sin(theta_y), 0, cos(theta_y)}
    };
    float Rz[3][3] = {
        {cos(theta_z), -sin(theta_z), 0},
        {sin(theta_z), cos(theta_z), 0},
        {0, 0, 1}
    };
    // rotate camera and up vectors
    // rotate about X
    vec_mat_mult(camera, Rx);
    vec_mat_mult(up, Rx);
    // rotate about Y
    vec_mat_mult(camera, Ry);
    vec_mat_mult(up, Ry);
    // rotate about Z
    vec_mat_mult(camera, Rz);
    vec_mat_mult(up, Rz);
    // move & scale camera
    camera[0] = 1.5 * E * camera[0] + center[0];
    camera[1] = 1.5 * E * camera[1] + center[1];
    camera[2] = 1.5 * E * camera[2] + center[2];
    printf("up: x=%f, y=%f, z=%f\n", up[0], up[1], up[2]);
    printf("camera: x=%f, y=%f, z=%f\n", camera[0], camera[1], camera[2]);
}

bbox_coords_t *calc_3D_coords(float *camera, float *up, float *center, float E) {
    int i;
    float a;
    bbox_coords_t *c = (bbox_coords_t *)calloc(1, sizeof(bbox_coords_t));
    c->top = (float *)calloc(3, sizeof(float));
    c->bottom = (float *)calloc(3, sizeof(float));
    c->topleft = (float *)calloc(3, sizeof(float));
    float tmp[3] = {center[0] - camera[0], center[1] - camera[1], center[2] - camera[2]};
    c->left = cross_product(up, tmp); // returns allocated mem
    c->right = cross_product(tmp, up); // returns allocated mem
    a = sqrt(pow(c->left[0], 2) + pow(c->left[1], 2) + pow(c->left[2], 2));

    for (i = 0; i < 3; i++) {
        c->left[i] = E / (2 * a) * c->left[i] + center[i];
        c->right[i] = E / (2 * a) * c->right[i] + center[i];
        c->top[i] = E / 2 * up[i] + center[i];
        c->bottom[i] = -E / 2 * up[i] + center[i];
        c->topleft[i] = E / 2 * up[i] + c->left[i];
    }

    return c;
}

void draw_image(unsigned char *img, float *z_buf, int rows, int cols, float **vertices, int num_vertices, int **triangles, int num_tris, bbox_coords_t *cds, float *camera) {
    int i, j, r, c;
    float D, n, d, dot1, dot2, dot3;
    float *plane, *cp1, *cp2, plane_neg[3];
    float image[3], intersect[3], v0[3], v1[3], v2[3];
    float tmp1[3], tmp2[3], tmp3[3], tmp4[3];

    // for each pixel in the image
    for (r = 0; r < rows; r++) {
        for (c = 0; c < cols; c++) {
            // calculate image pixel vector coordinates
            for (j = 0; j < 3; j++) {
                image[j] = cds->topleft[j] + ((float)c / (cols - 1) * (cds->right[j] - cds->left[j])) +
                           ((float)r / (rows - 1) * (cds->bottom[j] - cds->top[j]));
            }

            // for each triangle in the file
            for (i = 0; i < num_tris; i++) {
                // find plane eqn containing traingle
                for (j = 0; j < 3; j++) {
                    v0[j] =  vertices[triangles[i][0]][j];
                    v1[j] =  vertices[triangles[i][1]][j];
                    v2[j] =  vertices[triangles[i][2]][j];
                    tmp1[j] = v1[j] - v0[j];
                    tmp2[j] = v2[j] - v0[j];
                }

                plane = cross_product(tmp1, tmp2); // allocates mem

                for (j = 0; j < 3; j++)
                    plane_neg[j] = -1 * plane[j];

                D = dot_product(plane_neg, v0, 3);

                // find distance along image pixel ray to triangle (n/d)
                n = dot_product(plane_neg, camera, 3) - D;

                for (j = 0; j < 3; j++)
                    tmp1[j] = image[j] - camera[j];

                d = dot_product(plane, tmp1, 3);
                free(plane);

                // skip if ray is parallel to triangle
                if (fabs(d) < 0.001)
                    continue;

                // find ray-plane intersection
                for (j = 0; j < 3; j++)
                    intersect[j] = camera[j] + (n / d) * (image[j] - camera[j]);

                // determine if intersection lies within triangle
                for (j = 0; j < 3; j++) {
                    tmp1[j] = v2[j] - v0[j];
                    tmp2[j] = v1[j] - v0[j];
                    tmp3[j] = intersect[j] - v0[j];
                    tmp4[j] = v1[j] - v0[j];
                }

                cp1 = cross_product(tmp1, tmp2); // allocates mem
                cp2 = cross_product(tmp3, tmp4); // allocates mem
                dot1 = dot_product(cp1, cp2, 3);
                free(cp1);
                free(cp2);

                for (j = 0; j < 3; j++) {
                    tmp1[j] = v0[j] - v1[j];
                    tmp2[j] = v2[j] - v1[j];
                    tmp3[j] = intersect[j] - v1[j];
                    tmp4[j] = v2[j] - v1[j];
                }

                cp1 = cross_product(tmp1, tmp2); // allocates mem
                cp2 = cross_product(tmp3, tmp4); // allocates mem
                dot2 = dot_product(cp1, cp2, 3);
                free(cp1);
                free(cp2);

                for (j = 0; j < 3; j++) {
                    tmp1[j] = v1[j] - v2[j];
                    tmp2[j] = v0[j] - v2[j];
                    tmp3[j] = intersect[j] - v2[j];
                    tmp4[j] = v0[j] - v2[j];
                }

                cp1 = cross_product(tmp1, tmp2); // allocates mem
                cp2 = cross_product(tmp3, tmp4); // allocates mem
                dot3 = dot_product(cp1, cp2, 3);
                free(cp1);
                free(cp2);

                if (dot1 < 0 || dot2 < 0 || dot3 < 0)
                    continue;

                // check distance to triangle
                if ((n / d) > z_buf[r * cols + c])
                    continue;
                else
                    z_buf[r * cols + c] = n / d;

                // set pixel color
                img[r * cols + c] = 155 + (i % 100);
            }
        }

        printf("\rRendering %d%%", (int)(r / 255.0 * 100));
        fflush(stdout);
    }

    printf("\n");
}

void write_image(unsigned char *img, int rows, int cols, char *filename) {
    FILE *fpt = fopen(filename, "w");
    fprintf(fpt, "P5 %d %d 255\n", cols, rows);
    fwrite(img, cols * rows, 1, fpt);
    fclose(fpt);
}

void vec_mat_mult(float *vector, float matrix[3][3]) {
    int i, j;
    float tmp[3] = {0};

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            tmp[i] += vector[j] * matrix[j][i];

    for (i = 0; i < 3; i++)
        vector[i] = tmp[i];
}

float dot_product(float *v1, float *v2, int n) {
    int i;
    float result = 0.0;

    for (i = 0; i < n; i++)
        result += v1[i] * v2[i];

    return result;
}

float *cross_product(float *v1, float *v2) {
    float *result = (float *)calloc(3, sizeof(float));
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return result;
}
