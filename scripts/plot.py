from matplotlib import pyplot as plt
import numpy as np

x = [83, 
68.5619, 
52.3987, 
35.1461, 
17.4397, 
-0.0848067, 
-16.7918, 
-32.0456, 
-45.2106, 
-55.651, 
-62.7313, 
-65.8157, 
-64.5009, 
-59.3121, 
-51.0069, 
-40.3429, 
-28.0776, 
-14.9685, 
-1.77318, 
10.7508, 
21.8458, 
30.7545, 
36.7191, 
38.9822, 
37.0164, 
31.2144, 
22.1992, 
10.5937, 
-2.97902, 
-17.8961, 
-33.5347, 
-49.2717, 
-64.4843, 
-78.5495, 
-90.8444, 
-100.746, 
-107.762, 
-111.921, 
-113.382, 
-112.306, 
-108.85, 
-103.176, 
-95.4418, 
-85.8071, 
-74.4314, 
-61.474, 
-47.0943, 
-31.4518, 
-14.7057, 
2.94069, 
21.1099, 
39.3806, 
57.3314, 
74.5411, 
90.5883, 
105.052, 
117.51, 
127.542, 
134.727, 
138.642, 
138.867, 
135.144, 
127.871, 
117.608, 
104.917, 
90.3594, 
74.4953, 
57.8862, 
41.0933, 
24.6249, 
8.77838, 
-6.20181, 
-20.0711, 
-32.585, 
-43.4989, 
-52.5682, 
-59.5486, 
-64.1953, 
-66.3612, 
-66.2884, 
-64.3161, 
-60.7837, 
-56.0305, 
-50.3959, 
-44.2193, 
-37.7863, 
-31.1682, 
-24.3828, 
-17.4478, 
-10.3809, 
-3.19968, 
4.07802, 
11.4345, 
18.8522, 
26.3132, 
33.8, 
41.2948, 
48.7799, 
56.2376, 
63.6502]

y = [ -5, 
26.1856, 
48.0541, 
61.8914, 
68.9833, 
70.6158, 
68.0748, 
62.646, 
55.6155, 
48.2689, 
41.8923, 
37.7715, 
36.8988, 
39.0929, 
43.8785, 
50.7808, 
59.3247, 
69.0352, 
79.4373, 
90.056, 
100.416, 
110.043, 
118.461, 
125.196, 
129.884, 
132.604, 
133.55, 
132.913, 
130.884, 
127.655, 
123.419, 
118.366, 
112.689, 
106.58, 
100.23, 
93.8303, 
87.5447, 
81.4204, 
75.4758, 
69.7296, 
64.2002, 
58.9063, 
53.8662, 
49.0986, 
44.622, 
40.4549, 
36.6159, 
33.1236, 
29.9963, 
27.2185, 
24.6374, 
22.0659, 
19.3171, 
16.2039, 
12.5392, 
8.13623, 
2.80779, 
-3.63309, 
-11.3734, 
-20.6002, 
-31.5005, 
-44.0823, 
-57.6378, 
-71.2801, 
-84.1223, 
-95.2777, 
-103.859, 
-108.98, 
-109.754, 
-105.583, 
-97.0327, 
-84.9561, 
-70.2077, 
-53.6416, 
-36.112, 
-18.4731, 
-1.57899, 
13.7161, 
26.7084, 
37.2948, 
45.5225, 
51.4389, 
55.0911, 
56.5264, 
55.7921, 
52.9794, 
48.3559, 
42.2333, 
34.9231, 
26.737, 
17.9866, 
8.98346, 
0.039225, 
-8.53448, 
-16.426, 
-23.3239, 
-28.9163, 
-32.8917, 
-34.9386, 
-34.7452]

z = [119, 
118.905, 
118.219, 
117.159, 
115.944, 
114.791, 
113.918, 
113.542, 
113.88, 
115.152, 
117.573, 
121.362, 
126.661, 
133.311, 
141.08, 
149.733, 
159.036, 
168.757, 
178.661, 
188.514, 
198.083, 
207.135, 
215.435, 
222.751, 
228.888, 
233.816, 
237.544, 
240.083, 
241.442, 
241.63, 
240.658, 
238.535, 
235.27, 
230.874, 
225.355, 
218.724, 
211.021, 
202.405, 
193.065, 
183.193, 
172.978, 
162.61, 
152.279, 
142.174, 
132.486, 
123.405, 
115.12, 
107.822, 
101.7, 
96.8878, 
93.2923, 
90.7641, 
89.1537, 
88.3114, 
88.0879, 
88.3336, 
88.8989, 
89.6344, 
90.3906, 
91.0179, 
91.3668, 
91.3371, 
91.0256, 
90.5786, 
90.1421, 
89.8624, 
89.8856, 
90.3578, 
91.4253, 
93.2061, 
95.7062, 
98.9037, 
102.776, 
107.302, 
112.459, 
118.225, 
124.579, 
131.497, 
138.922, 
146.651, 
154.444, 
162.063, 
169.269, 
175.821, 
181.482, 
186.07, 
189.633, 
192.28, 
194.117, 
195.252, 
195.791, 
195.842, 
195.511, 
194.907, 
194.136, 
193.304, 
192.52, 
191.891, 
191.523, 
191.524
]

# mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(x, y, z, label='splined curve')
ax.legend()
plt.show()