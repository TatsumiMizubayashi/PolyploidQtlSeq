namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ作成設定値 インターフェース
    /// </summary>
    public interface IGraphSettingValue
    {
        /// <summary>
        /// グラフ画像の幅(Pixel)を取得する。
        /// </summary>
        int FigureWidth { get; }

        /// <summary>
        /// グラフ画像の高さ(Pixel)を取得する。
        /// </summary>
        int FigureHeight { get; }

        /// <summary>
        /// X軸の目盛り間隔(MB)を取得する。
        /// </summary>
        int XAxisMajorStep { get; }
    }
}
