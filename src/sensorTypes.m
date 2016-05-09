function st = armgenSensorTypes()

st = [];
st.inertial = 1; % (9D classical)
st.joint = 2; % direct joint value (1D)
st.globalpos = 3; % absolute position (3D)
st.globalz = 4; % global z axis vector (3D)
st.globalposAndz = 5; % special position + z axis (6D)