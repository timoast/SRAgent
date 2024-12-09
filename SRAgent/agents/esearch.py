# import
## batteries
import os
from typing import Annotated, List, Dict, Any, Callable
## 3rd party
from Bio import Entrez
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import HumanMessage, AIMessage
## package
from SRAgent.tools.esearch import esearch


# functions
def create_esearch_agent(model_name: str="gpt-4o-mini") -> Callable:
    """
    Create an agent that uses Entrez esearch to help complete a task.
    """
    model = ChatOpenAI(model_name=model_name, temperature=0.0)
    agent = create_react_agent(
        model=model,
        tools=[esearch],
        state_modifier="\n".join([
            "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
            "Based on the task provided by your supervisor, use Entrez esearch to help complete the task.",
            "If the sra or gds database does not return findings, try the other database.",
            "Provide a concise summary of your findings; use lists when possible; do not include helpful wording.",
        ])
    )
    
    @tool
    def invoke_esearch_agent(
        message: Annotated[str, "Message to the esearch agent"]
    ) -> Annotated[str, "Response from the esearch agent"]:
        """
        Invoke the esearch agent to perform a task.
        """
        # Invoke the agent with the message
        result = agent.invoke({"messages": [HumanMessage(content=message)]})
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="esearch_agent")]
        }
    return invoke_esearch_agent

if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    # test agent
    agent = create_esearch_agent()
    input = {"message": "I am looking for information about Entrez ID 35966237"}
    print(agent.invoke(input))

