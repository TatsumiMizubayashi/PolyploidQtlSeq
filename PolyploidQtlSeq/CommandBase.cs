using McMaster.Extensions.CommandLineUtils;

namespace PolyploidQtlSeq
{
    /// <summary>
    /// コマンドの基底クラス
    /// </summary>
    [HelpOption("-h|--help", Description = "Show help message.")]
    internal abstract class CommandBase
    {
        /// <summary>
        /// コマンドを実行する。
        /// </summary>
        /// <param name="app">App</param>
        /// <returns>終了コード</returns>
        public abstract Task<int> OnExecuteAsync(CommandLineApplication app);
    }
}
