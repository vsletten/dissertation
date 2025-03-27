# CLAUDE.md - Guidelines for Agentic Coding Assistants

## Build, Lint and Test Commands
- Build: `npm run build`
- Lint: `npm run lint`
- Test all: `npm run test`
- Test single: `npm test -- -t "test name"`
- Dev server: `npm run dev`

## Code Style Guidelines
- **Formatting**: Use Prettier with default settings
- **Imports**: Group imports (1. React/framework 2. External libraries 3. Internal modules)
- **Types**: Prefer TypeScript interfaces over types, explicit return types on functions
- **Naming**: camelCase for variables/functions, PascalCase for classes/components
- **Components**: Functional components with React hooks
- **Error Handling**: Use try/catch with specific error types and logging
- **Documentation**: JSDoc for public APIs, inline comments for complex logic
- **State Management**: Use React Context or dedicated state library for global state

Update this file as the project evolves.