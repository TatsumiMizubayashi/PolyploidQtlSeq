using Cysharp.Diagnostics;
using PolyploidQtlSeqCore.IO;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEff
    /// </summary>
    internal class SnpEff
    {
        private readonly SnpEffSettings _option;

        /// <summary>
        /// SnpEffを作成する。
        /// </summary>
        /// <param name="option">オプション</param>
        public SnpEff(SnpEffSettings option)
        {
            _option = option;
        }

        /// <summary>
        /// SnpEffを実行する。
        /// </summary>
        /// <param name="inputVcf">入力VCFファイル</param>
        /// <param name="outputVcfFilePath">出力VCFファイルのPath</param>
        /// <returns>SnpEff済みVCFファイル(圧縮されていない）</returns>
        public async ValueTask<VcfFile> RunAsync(VcfFile inputVcf, string outputVcfFilePath)
        {
            if (!_option.CanSneEff) throw new ArgumentException("Unable to execute SnpEff.");

            try
            {
                await RunSnpEffAsync(inputVcf, outputVcfFilePath);

                return new VcfFile(outputVcfFilePath);
            }
            catch
            {
                throw;
            }
        }

        /// <summary>
        /// snpEffコマンドを作成する。
        /// </summary>
        /// <param name="inputVcf">入力VCF</param>
        /// <returns></returns>
        private string CreateCommand(VcfFile inputVcf)
        {
            var command = $"snpEff -Xms2g -Xmx{_option.MaxHeap.Value}g ";
            if (_option.ConfigFile.HasFile) command += $"-c {_option.ConfigFile.Path} ";
            command += $"{_option.Database.Value} {inputVcf.Path} -noStats";

            return command;
        }

        /// <summary>
        /// SnpEffを実行する。
        /// </summary>
        /// <param name="inputVcf">入力VCF</param>
        /// <param name="outputVcfFilePath">出力VCFファイルPath</param>
        /// <returns></returns>
        private async ValueTask RunSnpEffAsync(VcfFile inputVcf, string outputVcfFilePath)
        {
            var command = CreateCommand(inputVcf);
            CommandLog.Add(command);

            try
            {
                using var writer = new StreamWriter(outputVcfFilePath);
                await foreach (var line in ProcessX.StartAsync(command))
                {
                    writer.WriteLine(line);
                }
            }
            catch (ProcessErrorException ex)
            {
                Log.AddRange(ex.ErrorOutput);

                // await foreach()はExitCodeが0以外か標準エラーがある場合例外を発生する。
                // 処理はうまくいっているが警告がある場合標準エラーに追加され例外が発生してしまう。
                // ExitCodeが0ならば問題ないと判断し終了する。
                if (ex.ExitCode == 0) return;

                throw;
            }
        }
    }
}
