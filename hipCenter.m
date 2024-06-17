function  [right_hip_center, left_hip_center, center] = hipCenter(right,left,sacrum)

    v1 = norm(right-left);

    vector3 = (right - left)/v1;

    v2 = norm(left-sacrum);

    fake_left = (left - sacrum)/v2;

    v3 = norm(right-sacrum);

    fake_right = (right - sacrum)/v3;

    vector1 = cross(fake_left, fake_right);

    vector2 = cross(vector3, vector1);

    center = horzcat(vector1', vector2', vector3');

    
% find center
    ASIS_breadth = norm(right - left); % ASIS breadth

    origin = (right + left)/2;

    breadth_vec_right = [0.24*ASIS_breadth; -0.21*ASIS_breadth; 0.32*ASIS_breadth];

    new_vec_right = center*breadth_vec_right;

    right_hip_center = origin' + new_vec_right;

    breadth_vec_left = [0.24*ASIS_breadth; -0.21*ASIS_breadth; -0.32*ASIS_breadth];

    new_vec_left = center*breadth_vec_left;

    left_hip_center = origin' + new_vec_left;

end