namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// スライディングウインドウ解析設定
    /// </summary>
    internal class SlidingWindowAnalysisSettings
    {

        /// <summary>
        /// スライディングウインドウ解析設定を作成する。
        /// </summary>
        /// <param name="settingValue">設定値値</param>
        public SlidingWindowAnalysisSettings(ISlidingWindowAnalysisSettingValue settingValue)
        {
            WindowSize = new WindowSize(settingValue.WindowSize);
            StepSize = new StepSize(settingValue.StepSize);
        }

        /// <summary>
        /// Window Size
        /// </summary>
        public WindowSize WindowSize { get; }

        /// <summary>
        /// StepSize
        /// </summary>
        public StepSize StepSize { get; }
    }
}
