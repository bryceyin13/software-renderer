# A Software renderer in C ++

### Intro

A software renderer in C++.

### Features:

- MVP transforms
- Triangle rasterization & Back face culling
- Texture mapping
- Z-buffer
- Perspective projection
- Programmable vertex & fragment shader
- Implementation of SSAO

### Requirements

- C++ (>=14)

### Build & Run

It is implemented and tested on Windows with VS2022. 

- Use the following command to get source code:

```
> git clone https://github.com/bryceyin13/Software-Renderer
```

- Create a new empty project with Visual Studio and name it whatever you like.
- In Solution Explorer, add all .h files to *Header files* and all .cpp files to *Source files*.
- Build & run.

### Result

shadow mapping

![shadow](https://github.com/bryceyin13/Software-Renderer/blob/main/image/shadow.png)

SSAO + shadow mapping (See the shallow shadow around the man's eyes, nose and mouth.)

![SSAO](https://github.com/bryceyin13/Software-Renderer/blob/main/image/shadow%2BSSAO.png)

### References

- [PBRT](https://pbr-book.org/4ed/contents)
- [GAMES online courses](https://games-cn.org/gamescoursescollection/)
- Ssloy's guide of [tinyrenderer](https://github.com/ssloy/tinyrenderer).