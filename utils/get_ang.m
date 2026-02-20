function ang = get_ang(vec1, vec2)
    tmp1 = vec1 / ( norm(vec1) + 1e-16 );
    tmp2 = vec2 / ( norm(vec2) + 1e-16 );
    cos_ang = tmp1' * tmp2;
    cos_ang = min(cos_ang, 1.0 - 1e-16);
    cos_ang = max(cos_ang, -1.0 + 1e-16);
    ang = acos(cos_ang);
end