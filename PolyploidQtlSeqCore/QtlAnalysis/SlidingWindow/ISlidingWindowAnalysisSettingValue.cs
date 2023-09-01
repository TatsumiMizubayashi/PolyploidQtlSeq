namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// スライディングウインドウ解析設定値インターフェイス
    /// </summary>
    public interface ISlidingWindowAnalysisSettingValue
    {
        /// <summary>
        /// Windowサイズ(kbp)を取得する。
        /// </summary>
        int WindowSize { get; }

        /// <summary>
        /// Window移動量(kbp)を取得する。
        /// </summary>
        int StepSize { get; }
    }
}
