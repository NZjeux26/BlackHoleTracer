# Blackhole Raytraccing Simulator

This project started off as a simple attempt to replicate the Lumiet Blackhole image from his 1978 Paper. Instead of using complicated shaders and C++ I wanted to use just SDL and C to replicate the image, and since I wanted this to jsut be a 2D image and not a full 3D simulation I thought it would be much simpler achievable even without LLM help.

It wasn't and now I have a 3D simulation of a Blackhole in OpenGL with GLSL and C as the backing langauage.

LLM assitance was used extensively, I have no real calculus background and I was using tools (after leaving SDL) that were beyond my skills to learn so rapidly with such complicated physics. Claude was the main LLM used, with assistance with the physics side from chatGPT (Graduate Physics trained submodel) and research help with Deep Research from Gemini. 

This however didn't mean my own programmer instinct wasn't used, LLMs work best with very directed questions and thus I did all code structutre and project layout along with the hours of debugging trying to find out why certian things happen as well as actual coding where needed to chnage things, bnut this was more in the style of working on an already established codebase vs blasting out code from scratch. 

I also worked at a VFX stuido, so leveraging the very talented crewmates I have to offer insight in either the physics, lighting, tools/language and just general ideas to either fix things on my own, or formulate better questions for the LLMs.

This is the current up to date production. At the moment it's just a schwarzschild blackhole but in te process of making into a Kerr blackhole. See the Images folder for more renders in different skyboxes, camera angles and directions.

![blackhole_gpu](https://github.com/user-attachments/assets/a5491207-b468-4649-9b57-963138ae6bd0)
