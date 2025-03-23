# import
## batteries
import os
import asyncio
from typing import Annotated, List, Optional, Callable
## 3rd party
from pydantic import BaseModel, Field
from Bio import Entrez
from langchain_core.tools import tool
from langgraph.prebuilt import create_react_agent
from langchain_core.runnables.config import RunnableConfig
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
## package
from SRAgent.agents.utils import set_model
from SRAgent.agents.tissue_ontology import create_tissue_ontology_agent

# classes
class UBERON_IDS(BaseModel):
    #descriptions: List[str] = Field(description="Each unique tissue description provided in the input.")
    ids: List[str] = Field(description="The selected Uberon ID (UBERON:XXXXXXX)")

# functions
def create_tissue_ontology_workflow(
    model_name: Optional[str]=None,
    return_tool: bool=True,
) -> Callable:
    # create model
    model = set_model(model_name=model_name, agent_name="tissue_ontology")

    # set tools
    tools = [
        create_tissue_ontology_agent(),
    ]
  
    # state modifier
    state_mod = "\n".join([
        "# Introduction",
        " - You are a helpful senior bioinformatician assisting a researcher with a task involving classifying one or more tissue.",
        " - You will be provided with a free text description of the tissues.",
        " - Your task is to categorize the tissues based on the Uberon ontology.",
        " - You must find the single most suitable Uberon ontology term that best describes the tissue description.",
        " - You have a set of tools that can help you with this task.",
        "# Tools",
        " - create_tissue_ontology_agent: Use this tool to find the most suitable Uberon ontology term that best describes the tissue description.",
        "# Workflow",
        " 1. Identify each unique tissue description in the input.",
        "   - For example, 'brain cortex; eye lens; aortic valve;' should be split into the following separate descriptions:",
        "     - brain cortex",
        "     - eye lens",
        "     - aortic valve",
        " 2. For each description (e.g., \"brain cortex\"), use the create_tissue_ontology_agent tool to find the most suitable Uberon ontology term.",
        "   - You MUST use the create_tissue_ontology_agent tool for EACH tissue description."
        "# Response",
        " - Provide the list of most suitable Uberon ontology IDs (UBERON:XXXXXXX) that best describe each tissue description.",
    ])
    # create agent
    agent = create_react_agent(
        model=model,
        tools=tools,
        state_modifier=state_mod,
        response_format=UBERON_IDS,
    )

    # create tool
    @tool
    async def invoke_tissue_ontology_workflow(
        messages: Annotated[List[BaseMessage], "Messages to send to the Tissue Ontology agent"],
        config: RunnableConfig,
    ) -> Annotated[dict, "Response from the Tissue Ontology agent"]:
        """
        Invoke the Tissue Ontology agent with a message.
        The Tissue Ontology agent will annotate each tissue description with the most suitable Uberon term.
        """
        # The React agent expects messages in this format
        response = await agent.ainvoke({"messages" : messages}, config=config)
        # return ontology term IDs
        return response['structured_response'].ids
    return invoke_tissue_ontology_workflow

# main 
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv(override=True)

    async def main():
        # create workflow
        workflow = create_tissue_ontology_workflow()

        # Example 1: Complex tissue description example
        print("\n=== Example 1: Complex tissue description example ===")
        msg = "Categorize the following tissues: the thin layer of epithelial cells lining the alveoli in lungs; brain cortex; eye lens"
        input = {"messages": [HumanMessage(content=msg)]}
        results = await workflow.ainvoke(input)
        print(results)

    # run
    asyncio.run(main())