import sys
import asyncio
from typing import List, Dict, Any
from langchain_openai import ChatOpenAI
from langchain_core.prompts import PromptTemplate

def set_model(model_name: str="gpt-4o-mini", temperature: float=0.1, reasoning_effort: str="low") -> ChatOpenAI:
    if  model_name.startswith("gpt-4o"):
        model = ChatOpenAI(model_name=model_name, temperature=temperature, reasoning_effort=None)
    elif model_name.startswith("o3") or model_name.startswith("o1"):
        model = ChatOpenAI(model_name=model_name, temperature=None, reasoning_effort=reasoning_effort)
    else:
        raise ValueError(f"Model {model_name} not supported")
    return model

def create_step_summary_chain(model: str="gpt-4o-mini", max_tokens: int=45):
    """
    Create a chain of tools to summarize each step in a workflow.
    Args:
        model: The OpenAI model to use for the language model.
        max_tokens: The maximum number of tokens to use for the summary.
    Returns:
        A chain of tools to summarize each step in a workflow.
    """
    # Create the prompt template
    template = "\n".join([
        "Concisely summarize the provided step in the langgraph workflow.",
        f"The summary must be {max_tokens} tokens or less.",
        "Do not use introductory words such as \"The workflow step involves\"",
        "Write your output as plain text instead of markdown.",
        "#-- The workflow step --#",
        "{step}"
    ])
    prompt = PromptTemplate(input_variables=["step"], template=template)

    # Initialize the language model
    llm = ChatOpenAI(model_name=model, temperature=0, max_tokens=max_tokens)

    # Return the LLM chain
    return prompt | llm


async def create_agent_stream(
    input,  
    create_agent_func,
    config: dict={}, 
    summarize_steps: bool=False
) -> str:
    """
    Create an Entrez agent and stream the steps.
    Args:
        input: Input message to the agent.
        create_agent_func: Function to create the agent.
        config: Configuration for the agent.
        summarize_steps: Whether to summarize the steps.
    Returns:
        The final step message.
    """
    # create entrez agent
    agent = create_agent_func(return_tool=False)

    # create step summary chain
    step_summary_chain = create_step_summary_chain() if summarize_steps else None
    
    # invoke agent
    step_cnt = 0
    final_step = ""
    async for step in agent.astream(input, stream_mode="values", config=config):
        step_cnt += 1
        final_step = step
        # summarize step
        if step_summary_chain:
            msg = step_summary_chain.invoke({"step": step.get("messages")})
            print(f"Step {step_cnt}: {msg.content}", file=sys.stderr)
        else:
            try:
                if "messages" in step and step["messages"]:
                    last_msg = step["messages"][-1].content
                    if last_msg != "":
                        print(f"Step {step_cnt}: {last_msg}", file=sys.stderr)
                    else:
                        step_cnt -= 1
            except (KeyError, IndexError, AttributeError):
                print(f"Step {step_cnt}: {step}", file=sys.stderr)
    try:
        final_step = final_step["agent"]["messages"][-1].content
    except KeyError:
        try:
            final_step = final_step["messages"][-1].content
        except (KeyError, IndexError, AttributeError):
            if isinstance(final_step, str):
                return final_step
            return str(final_step)
    return final_step